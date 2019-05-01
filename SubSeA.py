#!/usr/bin/env python3                                          
import sys
import os

from collections import defaultdict
from operator import itemgetter
import itertools

import subprocess
import numpy as np

import scipy.stats
import scipy.special
import requests
import json
BASE_PATH='/scratch/asl47/PDB/'

def needleAlign(pdb_1,chain_1,pdb_2,chain_2,needle_EXEC='./needle'):
    for pdb, chain in ((pdb_1,chain_1),(pdb_2,chain_2)):
        try:
            subprocess.run('grep -i {0}_{1} -A1 all_fasta.txt > {2}FASTA/{0}_{1}.fasta.txt'.format(pdb,chain,BASE_PATH),shell=True,check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            raw_FASTA=requests.get('http://www.ebi.ac.uk/pdbe/entry/pdb/{}/fasta'.format(pdb)).text
            lines=raw_FASTA.split('\n')
            for idx in range(0,len(lines),2):
                if chain in lines[idx][lines[idx].rfind('|')+1:]:
                    with open(BASE_PATH+'FASTA/{}_{}.fasta.txt'.format(pdb,chain),'w') as file_:
                        file_.write('>{}_{}\n'.format(pdb,chain)+lines[idx+1])
                    break
            #print('need to wget and append, then recall')
    subprocess.run('{EXEC} {BP}FASTA/{0}_{1}.fasta.txt {BP}FASTA/{2}_{3}.fasta.txt -gapopen 10.0 -gapextend 0.5 -outfile {BP}NEEDLE/{0}_{1}_{2}_{3}.needle'.format(pdb_1,chain_1,pdb_2,chain_2,EXEC=needle_EXEC,BP=BASE_PATH),shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def makeTypes(pdb):
    similaritythresh = 2.0
    type_ = {}
    with open(BASE_PATH+'INT/{}.int'.format(pdb)) as file_:
        for line in file_:
            l = line.strip().split('\t')
            chain = l[0]
            res = l[1]
            if len(l) == 4:
                ll = l[3].split()
                if float(ll[2]) > 0:
                    type_[chain+res] = tuple(reversed(ll[0].split('_')))
            if len(l) > 4:
                jjj = []
                for j in l[3:]:
                    jj = j.split()
                    if float(jj[2]) > 0:
                        jjj.append((jj,float(jj[2])))
                if jjj != []:
                    if not(len(jjj) == 2 and abs(jjj[0][1]-jjj[1][1]) < similaritythresh and jjj[0][0][0] == jjj[1][0][3] and jjj[0][0][3] == jjj[1][0][0]): #SYMMETRIC INTERFACE: ALWAYS PICK FIRST COLUMN, OTHERWISE SORT
                        jjj.sort(key=itemgetter(1),reverse=True)
                    type_[chain+res] = tuple(reversed(jjj[0][0][0].split('_')))
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
                seqs[read]+=ll[2]
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
    for i in range(len(seq)):
        if seq[i] != '-':
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

def writePmatrix(pdbs,inter1a,inter1b,inter2a,inter2b):
    mat = defaultdict(int)
    key1, key2 = set(), set()
    
    for a,b,c,d in zip(inter1a,inter1b,inter2a,inter2b):
        mat[(a+b,c+d)] += 1
        key1.add(a+b)
        key2.add(c+d)

    with open(BASE_PATH+'PALIGN/{}_{}_{}_{}.pmatrix'.format(*pdbs),'w') as pmatrix_out:
        pmatrix_out.write('\t')

        ##write column headers
        pmatrix_out.write('\t'.join(map(str,sorted(key2)))+'\n')

        for i in sorted(key1):
            ##write row header
            pmatrix_out.write('{}\t'.format(i))
            ##write row values
            pmatrix_out.write('\t'.join(map(str,(mat[(i,j)] for j in sorted(key2))))+'\n')

        pmatrix_out.write('\n')
    
                
def generateAssistiveFiles(pdbs):
    needleAlign(*pdbs,needle_EXEC='/rscratch/asl47/needle')
    
    type1=makeTypes(pdbs[0])
    type2=makeTypes(pdbs[2])
    (seq1,seq2),align=readNeedle(pdbs)

    int_1A,int_1B=makeInteractions(type1,seq1,pdbs[1])
    int_2A,int_2B=makeInteractions(type2,seq2,pdbs[3])



    writePialign(pdbs,int_1A,int_1B,int_2A,int_2B,align,seq1,seq2)
    writePmatrix(pdbs,int_1A,int_1B,int_2A,int_2B)
    

#if __name__ == "__main__":
#    if len(sys.argv) == 5:       
#        #needleAlign(*sys.argv[1:])
#        generateAssistiveFiles(sys.argv[1:])
#    else:
#        print('Wrong arguments')

def pcombine(pvalues):
    if 0 in pvalues:
        print(pvalues)
        sys.exit()
    return scipy.stats.combine_pvalues(pvalues,method='fisher')[1]

def binomialcdf(n,m1,m2,n1,n2,novg):
    return 1.0-scipy.stats.binom.cdf(n-1,novg,1.0*m1*m2/(n1*n2))

def MatAlign(pdb_1,chain_1,pdb_2,chain_2):
    matchmatrix = {}
    matchmatrixnorm = {}
    matchmatrixneedle = {}
    length = {}
    jkeys = set()
    kkeys = set()
    with open(BASE_PATH+'NEEDLE/{}_{}_{}_{}.needle'.format(pdb_1,chain_1,pdb_2,chain_2)) as needle_file:
        for line in needle_file:
            lss=line.strip().split()
            if 'Similarity' in line:
                if chain_1 not in matchmatrixneedle:
                    matchmatrixneedle[chain_1] = {}
                matchmatrixneedle[chain_1][chain_2] = line[line.find('(')+1:line.find('%')]
            elif 'Gaps' in line:
                gapped=line.strip().split()[2].split('/')
                noverlap=int(gapped[1])-int(gapped[0])
            elif line[0] != '#' and len(lss) == 4 and '_' in line:
                length[lss[0].upper()] = int(lss[-1])

    def readPmatrixFile():
        tsv = np.loadtxt(BASE_PATH+'PALIGN/{}_{}_{}_{}.pmatrix'.format(pdb_1,chain_1,pdb_2,chain_2),dtype=str,delimiter='\t')
        column_headers = tsv[0, 1:]
        row_headers = tsv[1:, 0]
        matrix = tsv[1:, 1:].astype(int)
        
        for swap in (column_headers,row_headers):
            if swap[0] == '  ':
                swap[0]='-'
        return column_headers, row_headers, matrix
                       

    column_headers, row_headers, matrix = readPmatrixFile()

    total = np.sum(matrix)
    rowsum = np.sum(matrix,axis=1)
    colsum = np.sum(matrix,axis=0)


    pvalues = []
    pmatrix = {}

    for row in range(len(matrix)):
        if row not in pmatrix:
            pmatrix[row] = {}
        for column in range(len(matrix[row])):
            if matrix[row,column] == 0 or row == 0 or column == 0:
                pmatrix[row][column] = 1.0
            else:
                pmatrix[row][column] = binomialcdf(matrix[row,column],rowsum[row],colsum[column],length['{}_{}'.format(pdb_1,chain_1)],length['{}_{}'.format(pdb_2,chain_2)],noverlap)

            if pmatrix[row][column] == 0:
                print('unknown')
                sys.exit()
            pvalues.append(pmatrix[row][column])
    pvalues.sort()


    ranked = defaultdict(list)
    for row in range(len(matrix)):
        for column in range(len(matrix[row])):
            ranked[pvalues.index(pmatrix[row][column])].append((row,column)) 

    already0 = set()
    already1 = set()
    matchmatrixnormtmp = []
    #print(ranked)
    for kk in sorted(ranked.keys()):
        for kkk in ranked[kk]:
            if kkk[0] not in already0 and kkk[1] not in already1:
                #fff.write('{}\t{}\t{}\t{}'.format(pbd_1,chain_1,pdb_2,chain_2)+'\t'+'\t'.join(list(h2[kkk[0]]))+'\t'+'\t'.join(list(h1[kkk[1]]))+'\t'+str(pmatrix[kkk[0]][kkk[1]])+'\t'+str(1.0*matrix[kkk[0]][kkk[1]]/total)+'\n')

                if row_headers[kkk[0]] != '-' and column_headers[kkk[1]] != '-':
                    if chain_1 not in matchmatrix:
                        matchmatrix[chain_1] = defaultdict(float)
                    if chain_1 not in matchmatrixnorm:
                        matchmatrixnorm[chain_1] = {}
                    matchmatrix[chain_1][chain_2] += matrix[kkk[0]][kkk[1]]/total
                    matchmatrixnormtmp.append(pmatrix[kkk[0]][kkk[1]])
                    jkeys.add(chain_1) # NECESSARY DUE TO GAPS
                    kkeys.add(chain_2) # NECESSARY DUE TO GAPS
                already0.add(kkk[0])
                already1.add(kkk[1])
    if matchmatrixnormtmp != []: 
        if chain_1 not in matchmatrixnorm:
            matchmatrixnorm[j] = {}#defaultdict(float)
        matchmatrixnorm[chain_1][chain_2] = pcombine(matchmatrixnormtmp)

    ##mass printing
    ##end


    if len(matchmatrix.keys()) > 0:
        maxperm = []
        minp = 1E6
        slk, slj = sorted([list(kkeys),list(jkeys)],key=len)

        total = scipy.special.perm(len(slj),len(slk),exact=True)
        if total < 1000000:
            for kk in itertools.permutations(slj,len(slk)):
                mapvalues = []
                needlevalues = []
                for kkk in range(len(slk)):
                    mapvalues.append(matchmatrixnorm[kk[kkk]][slk[kkk]])
                    needlevalues.append(matchmatrixneedle[kk[kkk]][slk[kkk]])
                    
                pc = pcombine(mapvalues) 
                if pc < minp:
                    minp = pc
                    maxperm = []
                    for kkk in range(len(slk)):
                        maxperm.append((kk[kkk],slk[kkk],mapvalues[kkk],needlevalues[kkk]))
        #else:
        #    print(i,total)

        mvtmp = []
        for kk in maxperm:
            #ffffff.write(str(i[0])+'\t'+str(i[1])+'\t')
            #for kkk in kk:
            #    ffffff.write(str(kkk)+'\t')
            mvtmp.append(kk[-2])
            #ffffff.write('\n')
            sout = '{}_{}\t{}_{}'.format(pdb_1,chain_1,pdb_2,chain_2)

            for kkk in kk:
                sout += str(kkk)+'\t'
            mvtmp.append(kk[-2])
            sout += '\n'
            return kk[2]
            #print(sout)
        #if mvtmp != [] and max(mvtmp) > thresh:
        #    fd.write('"'+str(i[0])+'" -> "'+str(i[1])+'";\n')
        
    
def runAlignBatch(pdb_code_file):
    results={}
    with open(pdb_code_file) as file_in:
        data=json.load(file_in)
        for pairing in data.values():
            for het, hom in itertools.product(*pairing):
                args=[het[:4].upper(),het[5],hom[:4].upper(),hom[5]]
                generateAssistiveFiles(args)
                results['{}_{}_{}_{}'.format(*args)]=MatAlign(*args)
    with open('pre_lim3.json', 'w') as file_:
          file_.write(json.dumps(results))
                

if __name__ == "__main__":
    try:
        runAlignBatch(sys.argv[1])
    except:
        print('wrong args')
        
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

def plotComps(ls):
    for l in ls:
        pass
            
def loadD():
    ds=[]
    with open('pre_lim.json') as file_in:
        data=json.load(file_in)
    for k,v in data.items():
        if v is None:
            print(k)
            continue
        d={}
        d['het']=k[:4]
        d['hom']=k[7:11]
        d['p']=v
        ds.append(d)
    return pd.DataFrame(ds)
                
    
