#!/usr/bin/env python3
import sys
import os

from collections import defaultdict, Counter
from operator import itemgetter
from time import sleep
import itertools

import subprocess
import numpy as np

import scipy.stats
import scipy.special
import requests
import json


BASE_PATH='/scratch/asl47/PDB/'

def pullFASTA(pdb,chain):
    for _ in range(3):
        try:
            lines = requests.get(f'http://www.ebi.ac.uk/pdbe/entry/pdb/{pdb}/fasta',timeout=5).text.split('\n')
            break
        except requests.Timeout:
            print('Timeout from EBI, pausing temporarily again')
            sleep(5)
    else:
        return False

    for idx in range(0,len(lines),2):
        if chain in lines[idx][lines[idx].rfind('|')+1:]:
            with open(f'{BASE_PATH}FASTA/{pdb}_{chain}.fasta.txt','w') as file_:
                file_.write('>{}_{}\n'.format(pdb,chain)+lines[idx+1])
                return True
    return False


    
##run needle alignment on two pbd_chain inputs
def needleAlign(pdb_1,chain_1,pdb_2,chain_2,needle_EXEC='./needle'):
    if os.path.exists(f'{BASE_PATH}NEEDLE/{pdb_1}_{chain_1}_{pdb_2}_{chain_2}.needle'):
        return True

    for pdb, chain in ((pdb_1,chain_1),(pdb_2,chain_2)):
        file_name = f'{BASE_PATH}FASTA/{pdb}_{chain}.fasta.txt'
        if os.path.exists(file_name) and os.path.getsize(file_name) > 20:
            continue

        ##see if fasta info is in given file
        try:
            subprocess.run(f'grep -i {pdb}_{chain} -A1 {BASE_PATH}FASTA/cleaned2_all_fasta.txt > {BASE_PATH}FASTA/{pdb}_{chain}.fasta.txt',shell=True,check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            ##could not find FASTA data in all_fasta file, try downloading
            print(f'Pulling FASTA data for {pdb}_{chain}')
            if not pullFASTA(pdb,chain):
                raise ValueError(f'Chain does not seem to exist for {pdb}_{chain}')

    subprocess.run(f'{needle_EXEC} {BASE_PATH}FASTA/{pdb_1}_{chain_1}.fasta.txt {BASE_PATH}FASTA/{pdb_2}_{chain_2}.fasta.txt -gapopen 10.0 -gapextend 0.5 -outfile {BASE_PATH}NEEDLE/{pdb_1}_{chain_1}_{pdb_2}_{chain_2}.needle',shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

def makeTypes(pdb):
    similaritythresh = 2.0
    type_ = {}
    with open(BASE_PATH+f'INT/{pdb}.int') as file_:
        for line_raw in file_:
            line = line_raw.split('\t')
            chain,res = line[:2]

            interactions=[(a,d,float(c)) for (a,b,c,d) in (Q.split() for Q in line[3:]) if c!='0']
                
            if interactions:
                if len(interactions) != 2 or abs(interactions[0][2]-interactions[1][2]) >= similaritythresh or interactions[0][0] != interactions[1][1] or interactions[0][1] != interactions[1][0]: 
                    interactions.sort(key=itemgetter(2),reverse=True)
                    
                type_[chain+res] = tuple(reversed(interactions[0][0].split('_')))
    return type_
        

##scrape FASTA sequence for both chains, and the alignment interactions 
def readNeedle(pdbs):
    seqs, align = ['', ''], ''
    read = 0
    similarity, score = None, None
#'/rscratch/asl47/liv_a.needle' or
    with open(BASE_PATH+'NEEDLE/{}_{}_{}_{}.needle'.format(*pdbs)) as file_:
        for line in file_:
            ##can skip empty or commented lines in .needle file
            if line == '\n' or line[0] == '#':
                if 'Similarity' in line:
                    similarity = float(line[line.find('(')+1 : line.find('%')])
                elif 'Score' in line:
                    score = float(line.split()[-1])
                continue
            
            ll = line.strip().split()
            ##add alternating sequence lines
            if len(ll) > 0 and ('_' in ll[0] or ll[0] == 'SEQUENCE'):
                seqs[read] += ll[2]
                read = not read
            ##if between chain A and chain B, record the alignment interactions
            elif read == 1:
                align += line[21:71].replace(' ','_').strip()

    assert similarity is not None, 'Probable error in needle calculation'
    return seqs, align, similarity, score

def makeInteractions(type_, seq, chain):

    def numalpha(s):
        if int(s) > 9:
            s = str(10+(int(s)-10)%26) # WRAP ROUND NUMBERS LARGER THAN 35 BACK TO a, b, c...
            s = chr(97+int(s)-10)
        return s

    c = 0
    inter_A = ''
    inter_B = ''
    for element in seq:
        if element != '-':
            chain_key = chain+str(c)
            if chain_key in type_:
                inter_A += numalpha(type_[chain_key][0])
                inter_B += type_[chain_key][1]
            else:
                inter_A += ' '
                inter_B += ' '
            c += 1
        else:
            inter_A += ' '
            inter_B += ' '
    return inter_A, inter_B

def writePialign(pdbs,inter1a,inter1b,inter2a,inter2b,align,seq1,seq2):
    iter_chunk=50
    with open(BASE_PATH+'PALIGN/{}_{}_{}_{}.pialign'.format(*pdbs),'w') as pialign_out:
        idx, idx_A, idx_B = 0,0,0
        for idx in range(0,max(len(seq1),len(seq2)),iter_chunk):
            slice_ = slice(idx,idx+iter_chunk)
            idx_A+=iter_chunk-seq1[slice_].count('-')
            idx_B+=iter_chunk-seq2[slice_].count('-')

            a=align[idx:idx+50]
           
            pialign_out.write(f'{inter1b[slice_]}\n{inter1a[slice_]}\n{seq1[slice_]}   {idx_A}\n{align[slice_]}\n{seq2[slice_]}   {idx_B}\n{inter2a[slice_]}\n{inter2b[slice_]}\n\n')

def makePmatrix(pdbs,inter1a,inter1b,inter2a,inter2b,write=True):
    matrix = defaultdict(int)
    rows, columns = set(), set()
    
    for a,b,c,d in zip(inter1a,inter1b,inter2a,inter2b):
        k1='-' if a+b == '  ' else a+b
        k2='-' if c+d == '  ' else c+d
        matrix[(k1,k2)] += 1
        rows.add(k1)
        columns.add(k2)

    rows, columns = sorted(rows), sorted(columns)
    
    if not write:
        matrix_rows=[]
        for i in rows:
            matrix_rows.append([matrix[(i,j)] for j in columns])
        return (rows, columns, np.array(matrix_rows,dtype=int)) 

    with open(BASE_PATH+'PALIGN/{}_{}_{}_{}.pmatrix'.format(*pdbs),'w') as pmatrix_out:
        ##write column headers
        pmatrix_out.write('\t'+'\t'.join(map(str,columns))+'\n')

        ##write row header + values
        for i in rows:            
            pmatrix_out.write(f'{i}\t'+'\t'.join(map(str,(matrix[(i,j)] for j in columns)))+'\n')

def readPmatrixFile(pdb_1,chain_1,pdb_2,chain_2):
    tsv = np.loadtxt(BASE_PATH+'PALIGN/{}_{}_{}_{}.pmatrix'.format(pdb_1,chain_1,pdb_2,chain_2),dtype=str,delimiter='\t')
    column_headers = tsv[0, 1:]
    row_headers = tsv[1:, 0]
    matrix = tsv[1:, 1:].astype(int)
    
    for swap in (column_headers,row_headers):
        if swap[0] == '  ':
            swap[0]='-'
    return column_headers, row_headers, matrix

    
                
def generateAssistiveFiles(pdbs,write_intermediates=False):
    needleAlign(*pdbs,needle_EXEC='/rscratch/asl47/needle')

    type1 = makeTypes(pdbs[0])
    type2 = makeTypes(pdbs[2])
    (seq1,seq2),align, similarity, score = readNeedle(pdbs)
 

    int_1A,int_1B = makeInteractions(type1,seq1,pdbs[1])
    int_2A,int_2B = makeInteractions(type2,seq2,pdbs[3])

    if write_intermediates:
        writePialign(pdbs,int_1A,int_1B,int_2A,int_2B,align,seq1,seq2)
        
    len_seq1,len_seq2 = (sum(c.isalpha() for c in seq) for seq in (seq1,seq2))
    no_overlap = len_seq1 + len_seq2 - len(align)
    
    return ((similarity, no_overlap,len_seq1,len_seq2,score), makePmatrix(pdbs,int_1A,int_1B,int_2A,int_2B,write_intermediates))

    


def pcombine(pvalues,comb_method='fisher'):
    if not pvalues or not all(pvalues):
        return 1
    else:
        if comb_method == 'tippett':
            return np.min(pvalues)
        return scipy.stats.combine_pvalues(pvalues,method=comb_method)[1]

def binomialcdf(n,m1,m2,n1,n2,novg,N_FACTOR=True):
    return scipy.stats.binom.sf(n-N_FACTOR,novg,m1*m2/(n1*n2))

def MatAlign(pdb_1,chain_1,pdb_2,chain_2,needle_result=None,matrix_result=None):
    if needle_result:
        similarity, noverlap, *length, score  = needle_result
        needle_length = sum(length) - noverlap
    else:
        (seq1,seq2),align, similarity, score = readNeedle([pdb_1,chain_1,pdb_2,chain_2])
        length = (sum(c.isalpha() for c in seq) for seq in (seq1,seq2))
        noverlap = sum(length) - len(align)
        needle_length = len(align)
        
        
    if matrix_result:
        row_headers,column_headers, matrix = matrix_result

    else:
        column_headers,row_headers, matrix = readPmatrixFile(pdb_1,chain_1,pdb_2,chain_2)


    ##trivially no possible matches
    if any(i<2 for i in matrix.shape):
        ##REJECT CODE
        return (1,1,1,0,similarity,score,needle_length,noverlap)

    #pmatrix = np.ones(matrix.shape)

    colsums,rowsums = np.meshgrid(*(np.sum(matrix,axis=I) for I in (1,0)),indexing='ij')
    

    pmatrix = binomialcdf(matrix,rowsums,colsums,*length,novg=noverlap,N_FACTOR=False)

    val_matrix=pmatrix[1:,1:].copy()

    def findMinElements(p_array,mins=[]):
        if np.min(p_array) <= 1:
            row_i,col_i = divmod(np.argmin(p_array),p_array.shape[1])
            mins.append(p_array[row_i,col_i])
            p_array[row_i,:] = np.inf
            p_array[:,col_i] = np.inf
            return findMinElements(p_array,mins)
        else:
            return mins

    alignment_scores = findMinElements(val_matrix)
    return (pcombine(alignment_scores,'fisher'),pcombine(alignment_scores,'stouffer'),pcombine(alignment_scores,'tippett'),len(alignment_scores),similarity,score,needle_length,noverlap)



        
from multiprocessing import Pool,Manager

def calculatePvalue(pdb_combination):
    (het,hom) = pdb_combination
    args=(het[:4].upper(),het[5],hom[:4].upper(),hom[5])
    try:
        n_r, m_r = generateAssistiveFiles(args)
        return (args,MatAlign(*args,needle_result=n_r,matrix_result=m_r))
    except Exception as e:
        print(het,hom,'!!Error!!:\t',e)
        return (args, 'error')
              
def paralleliseAlignment(pdb_pairs):
    print('Parellelising alignment')
    results = Manager().list()
    with Pool() as pool:
        for progress, (key,p_value) in enumerate(pool.imap_unordered(calculatePvalue,pdb_pairs,chunksize=50)):
            try:
                if os.path.exists('/scratch/asl47/PDB/NEEDLE/{}_{}_{}_{}.needle'.format(*key)):
                    os.remove('/scratch/asl47/PDB/NEEDLE/{}_{}_{}_{}.needle'.format(*key))
            except FileNotFoundError:
                print('File removal error')

            #results['{}_{}_{}_{}'.format(*key)]=p_value
            if p_value != 'error':
                results.append((key,)+p_value)
                
            if progress and progress % 50000 == 0:
                print(f'done another 50k ({progress})')
    print('Finished parallel mapping')
    return list(results)
