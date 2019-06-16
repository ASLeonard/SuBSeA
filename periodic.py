import pandas
from domains import readDomains, invertDomains
from SubSeA import paralleliseAlignment, calculatePvalue
from numpy.random import randint
import sys


def readHeteromers(file_path='~/Downloads/PeriodicTable.csv'):
    return pandas.read_csv(file_path)

def domainStrictness():
    pass

def scrapePDBs(df):
    with open('period_pdb_codes.txt', 'w') as file_out:
          file_out.write(', '.join(df['PDB ID']))

def linkedProteinGenerator(df):
    key_interfaces = 'List of interface types (all identical subunits are given the same code)'
    key_reduced = 'Interfaces included in the reduced topology (1 - yes; 0 - no)'

    N_SAMPLE_LIMIT = 500
    inverted_homodimer_domains = invertDomains(readDomains('HCath'))
    ALL_domains = readDomains('CATH-B')
    chain_map = chainMap()

    result_rows = []

    #i=0
    i_count=0   
    for _, row in df.iterrows():
        pdb = row['PDB ID']

        #if pdb not in chain_map:
        #    print('Skipping',pdb)
        #    continue

        domains = ALL_domains[pdb]

        ##empty domain, no point continuing
        if not domains:
            continue
        
        interactions={edge for pair, active in zip(row[key_interfaces].split(','),row[key_reduced].split(',')) if active == '1' for edge in pair.split('-')}


        for i,pair in enumerate(interactions):
            #e1, e2 = pair.split('-')

            ##heteromeric interaction
            if True or e1 != e2:
                domain_overlap = 1 ##need to implement
                
                for edge in (pair,):#(e1, e2):


                    ##relabel domain
                    if pdb in chain_map and edge in chain_map[pdb]:
                        edge = chain_map[pdb][edge]

                    ##no information on this subunit for domains
                    if edge not in domains:
                        continue
                    #try:
                    #    edge = domain_edges[ordered_edges.index(edge)]
                    #except:
                    #    print('bad length on', pdb)
                    #    print(ordered_edges, domain_edges)
                    #    continue
                    if tuple(domains[edge]) in inverted_homodimer_domains:
                        for comp in inverted_homodimer_domains[tuple(domains[edge])]:
                            #i_count+=1
                            #if i_count > 1000:
                            #    return
                            yield (pdb + '_' + edge, '{}_{}'.format(*comp.split('_')))
                        
    return

def randomProteinSampler(df, N_SAMPLE_LIMIT=100000):

    key_interfaces = 'List of interface types (all identical subunits are given the same code)'
    key_reduced = 'Interfaces included in the reduced topology (1 - yes; 0 - no)'

    inverted_homodimer_domains = invertDomains(readDomains('HCath'))
    homodimer_set = [pdb for proteins in inverted_homodimer_domains.values() for pdb in proteins]
    heterodimer_set = []
    ALL_domains = readDomains('CATH-B')
    chain_map = chainMap()

    result_rows = []

    for _, row in df.iterrows():
        pdb = row['PDB ID']
        domains = ALL_domains[pdb]

        ##empty domain, no point continuing
        if not domains:
            continue
        
        interactions=[pair for pair, active in zip(row[key_interfaces].split(','),row[key_reduced].split(',')) if active == '1']

        for i,pair in enumerate(interactions):
            e1, e2 = pair.split('-')

            ##heteromeric interaction
            if e1 != e2:
                domain_overlap = 1 ##need to implement
                
                for edge in (e1, e2):

                    ##relabel domain
                    if pdb in chain_map and edge in chain_map[pdb]:
                        edge = chain_map[pdb][edge]

                    ##no information on this subunit for domains
                    if edge not in domains:
                        continue

                    if tuple(domains[edge]) in inverted_homodimer_domains:
                        heterodimer_set.append('{}_{}'.format(pdb,edge))

    
    for hetero_index, homo_index_ in zip(randint(0,len(heterodimer_set),N_SAMPLE_LIMIT),randint(0,len(homodimer_set),N_SAMPLE_LIMIT)):
        yield (heterodimer_set[hetero_index], homodimer_set[homo_index_])
    else:
        return

    

def fnc(d):
    (het,hom) = d
    args=(het[:4].upper(),het[5],hom[:4].upper(),hom[5])
    return '{}||{}'.format(args[0],args[2])
    
from multiprocessing import Pool
def qqq(gener):
    with Pool() as pool:
        results = pool.map(fnc,gener,50)
    
    return sum(1 for _ in results)

from collections import defaultdict
def chainMap():
    chainmap = defaultdict(dict)
    with open('chain_map.txt','r') as file_:
        for line in file_:
            (pdb, file_chain, pdb_chain) = line.rstrip().split('\t')
            chainmap[pdb][file_chain] = pdb_chain
    return dict(chainmap)

    fc = open('chain_map.txt')
    for line in fc:
        l = line.strip().split('\t')
        if l[0].upper() not in chainmap:
            chainmap[l[0].upper()] = {}
        chainmap[l[0].upper()][l[1]] = l[2]

    chainmapset = set(chainmap.keys())
    chainset = defaultdict(list)
    
    fc2 = open('pdb_chains')
    for line in fc2:
        l = line.strip().split('\t')
        chainset[l[0].upper()].append(set(l[1].split()))

    pfaml = {}
    fp = open('pdb_pfam_mapping.txt')
    fp.readline()
    for line in fp:
        l = line.strip().split('\t')
        if l[0] not in pfaml:
            pfaml[l[0]] = defaultdict(list)
        if l[0] in chainmap:
            for i in chainmap[l[0]]:
                for j in chainset[l[0]]:
                    if l[1] in j and chainmap[l[0]][i] in j:
                        pfaml[l[0]][i].append(l[4])
        else:
            pfaml[l[0]][l[1]].append(l[4])

    pfam = {}
    for i in pfaml:
        if len(pfaml[i]) > 0:
            pfam[i] = {}
            for j in pfaml[i]:
                pfam[i][j] = ';'.join(sorted(list(pfaml[i][j]))) #';'.join(set(list(pfaml[i][j])))

    return pfam

import json

def main():
    df =readHeteromers()
    gener = linkedProteinGenerator(df)
    #print("LEN",sum(1 for _ in gener))
    #gener = linkedProteinGenerator(df)
    dic = {}
    
    if len(sys.argv) == 1:
        print('Running sequentially')
        for i,pdb in enumerate(gener):
            #print(i)
            #print('On pdb ',pdb)
            #if i % 500 == 0:
            #    print('Now on',i)
            try:
                results = calculatePvalue(pdb)
            except Exception as err:
                print('Error on {}'.format(pdb),err)

            dic['{}_{}_{}_{}'.format(*results[0])] = results[1]
    else:
        dic = paralleliseAlignment(gener)

    with open('full_run_2.dict','w') as f_out:
        json.dump(dic,f_out)
        
if __name__ == "__main__":
    main()


