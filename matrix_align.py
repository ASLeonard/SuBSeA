#!/usr/bin/env python                                                           
import sys
import os
import random
import networkx
import operator
import matplotlib
import pylab
import math
from collections import defaultdict
import os.path
import itertools
import scipy
import scipy.stats # import binom,combine_pvalues

def pcombine(pvalues):
    if 0 in pvalues:
        print(pvalues)
        sys.exit()
    return scipy.stats.combine_pvalues(pvalues,method='fisher')[1]

def binomialcdf(n,m1,m2,n1,n2,novg):
    #print((n,m1,m2,n1,n2,novg))
    return 1.0-scipy.stats.binom.cdf(n-1,novg,1.0*m1*m2/(n1*n2))


cl = [[sys.argv[1],sys.argv[3],sys.argv[2],sys.argv[4]]]



fd = open('matrixextract_uncut_long_pisa.dot','w')
fd.write('digraph matrixextract\n{\n')
thresh = 0.5
ffffff = open('matrixextract_uncut_long_pisa.max','w')
fffff = open('matrixextract_uncut_long_pisa.matrix','w')
ffff = open('matrixextract_uncut_long_pisa.gaps','w')
fff = open('matrixextract_uncut_long_pisa.out','w')

def something(ID):
    matchmatrix = {}
    matchmatrixnorm = {}
    matchmatrixneedle = {}
    length = {}
    with open('NEEDLE/{}.needle'.format(ID)) as needle_file:
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
                
    pass

c = 0
for i in cl:

    j = i[2] #chain_1
    k = i[3] #chain_2
    #print(i,i[0] in chains,i[1] in chains)
    #if c%1000 == 0:
    #    print(c,pl)
    #c += 1
    if 1 == 1:#i[0] in chains and i[1] in chains:
        matchmatrix = {}
        matchmatrixnorm = {}
        matchmatrixneedle = {}
        jkeys = set([])
        kkeys = set([])
        if 1 == 1:#for j in chains[i[0]]:
            #print(j)
            if 1 == 1:#for k in chains[i[1]]:
                #print(k)
                fn = 'PALIGN/'+str(i[0])+'_'+str(j)+'_'+str(i[1])+'_'+str(k)+'.pmatrix'
                if os.path.isfile(fn):
                    length = {}
                    ffn = open('NEEDLE/'+fn[7:-8]+'.needle')
                    #print(fn[:-10])
                    for line in ffn:
                        lss = line.strip().replace('(','').split()
                        if 'Similarity' in line:
                            if j not in matchmatrixneedle:
                                matchmatrixneedle[j] = {}
                            matchmatrixneedle[j][k] = lss[-1][:-2].strip()

                        if 'Gaps' in line:
                            noverlap = int(lss[-2].split('/')[1])-int(lss[-2].split('/')[0])

                        if line[0] != '#' and len(lss) == 4 and '_' in line:
                            length[lss[0].split('_')[0].upper()+'_'+lss[0].split('_')[1]] = int(lss[-1])
                            #print(line)
                            print(lss[0].split('_')[0].upper()+'_'+lss[0].split('_')[1],int(lss[-1]))

                    ff = open(fn)
                    #print(fn)
                    h1 = ff.readline().strip().split()
                    firstline = ff.readline().strip().split()
                    ftmp = []
                    if len(firstline) > 0:
                        for kk in firstline[1-int(firstline[0].isdigit()):]:
                            ftmp.append(int(kk))
                        matrix = [ftmp[:]]
                        h2 = []
                        values = ftmp
                        for line in ff:
                            l = line.strip().split()
                            if len(l) > 0:
                                h2.append(l[0])
                                mtmp = []
                                #vtmp = []
                                for kk in l[1:]:
                                    mtmp.append(int(kk))
                                    values.append(int(kk))
                                matrix.append(mtmp)
                                #values.extend(vtmp)
                        if len(matrix) > len(h2):
                            h2.insert(0,'-')
                        if len(matrix[0]) > len(h1):
                            h1.insert(0,'-')
                        ff.close()
                        rowsum = {}
                        colsum = defaultdict(int)
                        for kk in range(len(matrix)):
                            rowsum[kk] = sum(matrix[kk])
                            for kkk in range(len(matrix[kk])):
                                colsum[kkk] += matrix[kk][kkk]
                        jvalues = []
                        jmatrix = {}
                        for kk in range(len(matrix)):
                            if kk not in jmatrix:
                                jmatrix[kk] = {}
                            for kkk in range(len(matrix[kk])):
                                jmatrix[kk][kkk] = 1.0*matrix[kk][kkk]/(rowsum[kk]+colsum[kkk]-matrix[kk][kkk])
                                jvalues.append(jmatrix[kk][kkk])
                        jvalues.sort(reverse=True)
                        pvalues = []
                        pmatrix = {}
                        for kk in range(len(matrix)):
                            if kk not in pmatrix:
                                pmatrix[kk] = {}
                            for kkk in range(len(matrix[kk])):
                                pmatrix[kk][kkk] = binomialcdf(matrix[kk][kkk],rowsum[kk],colsum[kkk],length[str(i[0])+'_'+str(j)],length[str(i[1])+'_'+str(k)],noverlap)
                                if matrix[kk][kkk] == 0 or kk == 0 or kkk == 0:
                                    pmatrix[kk][kkk] = 1.0
                                if pmatrix[kk][kkk] == 0:
                                    print('matrix,kk,kkk,matrix[kk][kkk],rowsum[kk],colsum[kkk],length[str(i[0])+str(j)],length[str(i[1])+str(k)],noverlap,binomialcdf',matrix,kk,kkk,matrix[kk][kkk],rowsum[kk],colsum[kkk],length[str(i[0])+'_'+str(j)],length[str(i[1])+'_'+str(k)],noverlap,binomialcdf(matrix[kk][kkk],rowsum[kk],colsum[kkk],length[str(i[0])+'_'+str(j)],length[str(i[1])+'_'+str(k)],noverlap))
                                    sys.exit()
                                pvalues.append(pmatrix[kk][kkk])
                        pvalues.sort()                    
                        total = sum(values) # NOT jvalues or pvalues - THIS IS CORRECT
                        ranked = defaultdict(list)
                        for kk in range(len(matrix)):
                            for kkk in range(len(matrix[kk])):
                                ranked[pvalues.index(pmatrix[kk][kkk])].append((kk,kkk)) #jvalues
                        already0 = set([])
                        already1 = set([])
                        matchmatrixnormtmp = []
                        for kk in sorted(ranked.keys()):
                            #print(kk)
                            for kkk in ranked[kk]:
                                #print(kkk)
                                if kkk[0] not in already0 and kkk[1] not in already1:# and jmatrix[kkk[0]][kkk[1]] > 0:
                                    fff.write(str(i[0])+'\t'+str(j)+'\t'+str(i[1])+'\t'+str(k)+'\t'+'\t'.join(list(h2[kkk[0]]))+'\t'+'\t'.join(list(h1[kkk[1]]))+'\t'+str(pmatrix[kkk[0]][kkk[1]])+'\t'+str(1.0*matrix[kkk[0]][kkk[1]]/total)+'\n')#'+str(1.0*matrix[kkk[0]][kkk[1]]/colsum[kkk[1]])+'\t'+str(1.0*matrix[kkk[0]][kkk[1]]/rowsum[kkk[0]])+'\t'+str(1.0*matrix[kkk[0]][kkk[1]]/total)+'\n')
                                    if h2[kkk[0]] != '-' and h1[kkk[1]] != '-':
                                        if j not in matchmatrix:
                                            matchmatrix[j] = defaultdict(float)
                                        if j not in matchmatrixnorm:
                                            matchmatrixnorm[j] = {}#defaultdict(float)
                                        matchmatrix[j][k] += 1.0*matrix[kkk[0]][kkk[1]]/total
                                        #if len(already0) == 0:
                                        matchmatrixnormtmp.append(pmatrix[kkk[0]][kkk[1]])
                                        #matchmatrixnorm[j][k] += (1.0*matrix[kkk[0]][kkk[1]]/total)*jmatrix[kkk[0]][kkk[1]]#(0.5*matrix[kkk[0]][kkk[1]]/colsum[kkk[1]]+0.5*matrix[kkk[0]][kkk[1]]/rowsum[kkk[0]])
                                        jkeys.add(j) # NECESSARY DUE TO GAPS
                                        kkeys.add(k) # NECESSARY DUE TO GAPS
                                    already0.add(kkk[0])
                                    already1.add(kkk[1])
                        #if j in matchmatrix and k in matchmatrix[j]:
                        #    matchmatrixnorm[j][k] /= matchmatrix[j][k]
                        if matchmatrixnormtmp != []: 
                            if j not in matchmatrixnorm:
                                matchmatrixnorm[j] = {}#defaultdict(float)
                            matchmatrixnorm[j][k] = pcombine(matchmatrixnormtmp)
                                    
                else: 
                    ffff.write(str(fn)+'\n')

        if len(matchmatrix.keys()) > 0:
            fffff.write('@\t'+str(i[0])+'\t'+str(i[1])+'\t'+str(len(jkeys))+'\t'+str(len(kkeys))+'\n\t')
            for kkk in sorted(kkeys):
                fffff.write(str(kkk)+'\t')
            fffff.write('\n')
            for kk in sorted(jkeys):
                fffff.write(str(kk)+'\t')
                for kkk in sorted(kkeys):
                    fffff.write(str(matchmatrix[kk][kkk])+'\t')
                fffff.write('\n')
            fffff.write('\n\t')
            for kkk in sorted(kkeys):
                fffff.write(str(kkk)+'\t')
            fffff.write('\n')
            for kk in sorted(jkeys):
                fffff.write(str(kk)+'\t')
                for kkk in sorted(kkeys):
                    fffff.write(str(matchmatrixnorm[kk][kkk])+'\t')
                fffff.write('\n')
            fffff.write('\n\t')
            for kkk in sorted(kkeys):
                fffff.write(str(kkk)+'\t')
            fffff.write('\n')
            for kk in sorted(jkeys):
                fffff.write(str(kk)+'\t')
                for kkk in sorted(kkeys):
                    fffff.write(str(matchmatrixneedle[kk][kkk])+'\t')
                fffff.write('\n')
            fffff.write('-------------------\n')

        # FIND MAX MATCH COMBINATION
        if len(matchmatrix.keys()) > 0:
            maxperm = []
            if len(jkeys) >= len(kkeys):
                maxperm = []
                #maxmap = -1
                minp = 1E6
                slj = sorted(list(jkeys))
                slk = sorted(list(kkeys))
                total = int(math.factorial(len(kkeys))*scipy.special.comb(len(slj),len(kkeys)))
                if total < 1000000:
                    #ccc = 0
                    for kk in itertools.permutations(slj,len(kkeys)):
                        #if ccc%1000 == 0:
                        #    print(ccc,total)
                        #ccc += 1
                        mapvalues = []
                        needlevalues = []
                        for kkk in range(len(kkeys)):
                            mapvalues.append(matchmatrixnorm[kk[kkk]][slk[kkk]])
                            needlevalues.append(matchmatrixneedle[kk[kkk]][slk[kkk]])
                        #if sum(mapvalues) > maxmap:
                            #maxmap = sum(mapvalues)
                        pc = pcombine(mapvalues) 
                        if pc < minp:
                            minp = pc
                            maxperm = []
                            for kkk in range(len(kkeys)):
                                maxperm.append((kk[kkk],slk[kkk],mapvalues[kkk],needlevalues[kkk]))
                else:
                    print(i,total)
            if len(jkeys) < len(kkeys):
                maxperm = []
                #maxmap = -1
                minp = 1E6
                slj = sorted(list(jkeys))
                slk = sorted(list(kkeys))
                total = int(math.factorial(len(jkeys))*scipy.special.comb(len(slk),len(jkeys)))
                if total < 1000000:
                    for kk in itertools.permutations(slk,len(jkeys)):
                        mapvalues = []
                        needlevalues = []
                        for kkk in range(len(jkeys)):
                            mapvalues.append(matchmatrixnorm[slj[kkk]][kk[kkk]])
                            needlevalues.append(matchmatrixneedle[slj[kkk]][kk[kkk]])
                        #if sum(mapvalues) > maxmap:
                            #maxmap = sum(mapvalues)
                        pc = pcombine(mapvalues) 
                        if pc < minp:
                            minp = pc
                            maxperm = []
                            for kkk in range(len(jkeys)):
                                maxperm.append((slj[kkk],kk[kkk],mapvalues[kkk],needlevalues[kkk]))
                else: 
                    print(i,total)

            mvtmp = []
            for kk in maxperm:
                ffffff.write(str(i[0])+'\t'+str(i[1])+'\t')
                for kkk in kk:
                    ffffff.write(str(kkk)+'\t')
                mvtmp.append(kk[-2])
                ffffff.write('\n')
                sout = ''
                sout += str(i[0])+'\t'+str(i[1])+'\t'
                for kkk in kk:
                    sout += str(kkk)+'\t'
                mvtmp.append(kk[-2])
                sout += '\n'
                print("Final result",sout)
            if mvtmp != [] and max(mvtmp) > thresh:
                fd.write('"'+str(i[0])+'" -> "'+str(i[1])+'";\n')

fff.close()
ffff.close()
fffff.close()
ffffff.close()
fd.write('}')
fd.close()
