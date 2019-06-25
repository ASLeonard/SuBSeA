#!/usr/bin/env python3
import sys
import os

from collections import defaultdict, Counter
from operator import itemgetter
import itertools

import subprocess
import numpy as np

import scipy.stats
import scipy.special
import requests
import json

BASE_PATH='/scratch/asl47/PDB/'

def pullFASTA(pdb,chain):
    lines = requests.get('http://www.ebi.ac.uk/pdbe/entry/pdb/{}/fasta'.format(pdb),timeout=1).text.split('\n')
    for idx in range(0,len(lines),2):
        if chain in lines[idx][lines[idx].rfind('|')+1:]:
            with open(BASE_PATH+'FASTA/{}_{}.fasta.txt'.format(pdb,chain),'w') as file_:
                file_.write('>{}_{}\n'.format(pdb,chain)+lines[idx+1])
                return True
    return False


    
##run needle alignment on two pbd_chain inputs
def needleAlign(pdb_1,chain_1,pdb_2,chain_2,needle_EXEC='./needle'):
    if os.path.exists('{BP}NEEDLE/{0}_{1}_{2}_{3}.needle'.format(pdb_1,chain_1,pdb_2,chain_2,BP=BASE_PATH)):
        return
    for pdb, chain in ((pdb_1,chain_1),(pdb_2,chain_2)):
        file_name = '{2}FASTA/{0}_{1}.fasta.txt'.format(pdb,chain,BASE_PATH)
        if os.path.exists(file_name) and os.path.getsize(file_name):
            continue

        ##see if fasta info is in given file
        try:
            subprocess.run('grep -i {0}_{1} -A1 {2}FASTA/clean_all_fasta.txt > {2}FASTA/{0}_{1}.fasta.txt'.format(pdb,chain,BASE_PATH),shell=True,check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print('Pulling FASTA data for {}_{}'.format(pdb,chain))
            pullFASTA(pdb,chain)

    subprocess.run('{EXEC} {BP}FASTA/{0}_{1}.fasta.txt {BP}FASTA/{2}_{3}.fasta.txt -gapopen 10.0 -gapextend 0.5 -outfile {BP}NEEDLE/{0}_{1}_{2}_{3}.needle'.format(pdb_1,chain_1,pdb_2,chain_2,EXEC=needle_EXEC,BP=BASE_PATH),shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

def makeTypes(pdb):
    similaritythresh = 2.0
    type_ = {}
    with open(BASE_PATH+'INT/{}.int'.format(pdb)) as file_:
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
    with open(BASE_PATH+'NEEDLE/{}_{}_{}_{}.needle'.format(*pdbs)) as file_:
        for line in file_:
            ##can skip empty or commented lines in .needle file
            if line == '\n' or line[0] == '#':
                continue
            
            ll = line.strip().split()
            ##add alternating sequence lines
            if len(ll) > 0 and ('_' in ll[0] or ll[0] == 'SEQUENCE'):
                seqs[read] += ll[2]
                read = not read
            ##if between chain A and chain B, record the alignment interactions
            elif read == 1:
                align += line[21:71].replace(' ','_').strip()

    return seqs, align
        

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
           
            pialign_out.write('{}\n{}\n{}   {}\n{}\n{}   {}\n{}\n{}\n\n'.format(inter1b[slice_],inter1a[slice_],seq1[slice_],idx_A,align[slice_],seq2[slice_],idx_B,inter2a[slice_],inter2b[slice_]))

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
            pmatrix_out.write('{}\t'.format(i)+'\t'.join(map(str,(matrix[(i,j)] for j in columns)))+'\n')

    
                
def generateAssistiveFiles(pdbs,write_intermediates=False):
    needleAlign(*pdbs,needle_EXEC='/rscratch/asl47/needle')

    type1=makeTypes(pdbs[0])
    type2=makeTypes(pdbs[2])
    (seq1,seq2),align=readNeedle(pdbs)

    int_1A,int_1B=makeInteractions(type1,seq1,pdbs[1])
    int_2A,int_2B=makeInteractions(type2,seq2,pdbs[3])

    if write_intermediates:
        writePialign(pdbs,int_1A,int_1B,int_2A,int_2B,align,seq1,seq2)
        
    len_seq1,len_seq2 = (sum(c.isalpha() for c in seq) for seq in (seq1,seq2))
    no_overlap = len_seq1 + len_seq2 - len(align)
    
    return ((len_seq1,len_seq2,no_overlap), makePmatrix(pdbs,int_1A,int_1B,int_2A,int_2B,write_intermediates))

    


def pcombine(pvalues):
    if not pvalues or not all(pvalues):
        return None
    else:
        return scipy.stats.combine_pvalues(pvalues,method='fisher')[1]

def binomialcdf(n,m1,m2,n1,n2,novg):
    return scipy.stats.binom.sf(n-1,novg,m1*m2/(n1*n2))

def MatAlign(pdb_1,chain_1,pdb_2,chain_2,needle_result=None,matrix_result=None):

    def readNeedleFile():
        with open(BASE_PATH+'NEEDLE/{}_{}_{}_{}.needle'.format(pdb_1,chain_1,pdb_2,chain_2)) as needle_file:
            length = {}
            for line in needle_file:
                lss=line.strip().split()
             
                if 'Gaps' in line:
                    gapped=line.strip().split()[2].split('/')
                    noverlap=int(gapped[1])-int(gapped[0])
                elif line[0] != '#' and len(lss) == 4 and '_' in line:
                    length[lss[0].upper()] = int(lss[-1])
        return tuple(length.values()), noverlap

    def readPmatrixFile():
        tsv = np.loadtxt(BASE_PATH+'PALIGN/{}_{}_{}_{}.pmatrix'.format(pdb_1,chain_1,pdb_2,chain_2),dtype=str,delimiter='\t')
        column_headers = tsv[0, 1:]
        row_headers = tsv[1:, 0]
        matrix = tsv[1:, 1:].astype(int)
        
        for swap in (column_headers,row_headers):
            if swap[0] == '  ':
                swap[0]='-'
        return column_headers, row_headers, matrix
                       
    if needle_result:
        length, noverlap = needle_result[:2], needle_result[2]
    else:
        length, noverlap = readNeedleFile()
        
    if matrix_result:
        row_headers,column_headers, matrix = matrix_result
    else:
        column_headers,row_headers, matrix = readPmatrixFile()

    ##trivially no possible matches
    if any(i<2 for i in matrix.shape):
        return None

        

    #total = np.sum(matrix)
    rowsum = np.sum(matrix,axis=1)
    colsum = np.sum(matrix,axis=0)

    
    pmatrix = np.empty(matrix.shape)

    for row in range(len(matrix)):
        for column in range(len(matrix[row])):
            if matrix[row,column] == 0 or row == 0 or column == 0:
                pmatrix[row,column] = 1.0
            else:
                pmatrix[row,column] = binomialcdf(matrix[row,column],rowsum[row],colsum[column],*length,novg=noverlap)


    val_matrix=pmatrix[1:,1:].copy()

    def findMinElements(tm,mins=[]):
        if len(mins)>1000:
            raise Exception('Recursed a lot of times, probably wrong')
        if np.min(tm) <= 1:
            ri,ci = divmod(np.argmin(tm),tm.shape[1])
            mins.append(tm[ri,ci])
            tm[ri,:]=2
            tm[:,ci]=2
            return findMinElements(tm,mins)
        else:
            return mins
            
        
    #return(pmatrix)
    return pcombine(findMinElements(val_matrix))
    #print(tuple((Counter(np.min(pmatrix[1:,1:],axis=0)) & Counter(np.min(pmatrix[1:,1:],axis=1))).elements()))
    ##find min pvalues for row/columns uniquely
    #return pcombine(tuple((Counter(np.min(pmatrix[1:,1:],axis=0)) & Counter(np.min(pmatrix[1:,1:],axis=1))).elements()))



        
from multiprocessing import Pool,Manager

def calculatePvalue(pdb_combination):
    (het,hom) = pdb_combination
    args=(het[:4].upper(),het[5],hom[:4].upper(),hom[5])
    try:
        n_r, m_r = generateAssistiveFiles(args)
        return (args,MatAlign(*args,needle_result=n_r,matrix_result=m_r))
    except Exception as e:
        print(het,hom,e)
        return (args, 'error')
    
def runParallelAlign(pdb_code_file):
    
    with open(pdb_code_file,'r') as file_in:
        data=json.load(file_in)
        
    results = Manager().dict()
        
    with Pool() as pool:
        for pairing in data.values():
            for (key,p_value) in pool.map(calculatePvalue,itertools.product(*pairing),chunksize=50):
                results['{}_{}_{}_{}'.format(*key)]=p_value

        
    with open('pre_lim4.json', 'w') as file_:
          file_.write(json.dumps(results.copy()))
          
def paralleliseAlignment(pdb_pairs):
    print('Parellelising alignment')
    results = Manager().dict()
    with Pool() as pool:
        for (key,p_value) in pool.imap_unordered(calculatePvalue,pdb_pairs,chunksize=50):
            #print(key)
            results['{}_{}_{}_{}'.format(*key)]=p_value

    return results.copy()
        
    

if __name__ == "__main__":
    #try:
    runParallelAlign(sys.argv[1])
    #except:
    #    print('wrong args')
        
import random
def randomArrange():
    data=json.load(open('domain_groups.json'))
    pdbs=([],[])
    for pairing in data.values():
        for i,pair in enumerate(pairing):
            pdbs[i].extend(pair)

    random_samples={}
    for i,pairing in enumerate(data.values()):
        random_samples[str(i)]=([],[])
        for j,pair in enumerate(pairing):
            for _ in range(len(pair)):
                random_samples[str(i)][j].append(pdbs[j].pop(random.randint(0,len(pdbs[j])-1)))
        
    return random_samples

import pandas as pd
import matplotlib.pyplot as plt

def plotComps(ls):
    plt.figure()
    labels={'':'Domains',2:'Random Shuffle',4:'Random Shuffle'}
    vx=[]
    for l in ls:
        df = loadD(l)
        ps = df['p']
        ns,_,_=plt.hist(ps,bins=100,density=True,histtype='step',label=labels[l])
        vx.append(ns)
    vx=np.array(vx)
    plt.plot(np.linspace(0,1,100),vx[0]/np.mean(vx[1:],axis=0),c='k',lw=2)
    plt.plot([0,1],[1,1],'k--')

    plt.yscale('log',nonposy='mask')
    plt.legend()
    plt.show(block=False)
            
def loadD(i):
    ds=[]
    with open('pre_lim{}.json'.format(i)) as file_in:
        data=json.load(file_in)
    for k,v in data.items():
        if v is None:
            continue
        d={}
        d['het']=k[:4]
        d['hom']=k[7:11]
        d['p']=v
        ds.append(d)
    return pd.DataFrame(ds)
                
    
