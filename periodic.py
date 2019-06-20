import pandas
from domains import readDomains, invertDomains, getUniqueHomodimerDomains
from SubSeA import paralleliseAlignment, calculatePvalue
from numpy.random import randint, choice
from collections import defaultdict, Counter
import sys
import json
import argparse
import matplotlib.pyplot as plt
import numpy as np

def readHeteromers(file_path='~/Downloads/PeriodicTable.csv'):
    return pandas.read_csv(file_path)

def domainStrictness():
    pass

def scrapePDBs(df):
    with open('period_pdb_codes.txt', 'w') as file_out:
          file_out.write(', '.join(df['PDB ID']))

def heteromerHomodimerOverlap(df,H_max=None,periodic=True):
    key_interfaces = 'List of interface types (all identical subunits are given the same code)'
    key_reduced = 'Interfaces included in the reduced topology (1 - yes; 0 - no)'

    inverted_homodimer_domains2 = invertDomains(readDomains('HCath'))
    unis=getUniqueHomodimerDomains(H_max)
    inverted_homodimer_domains = {k:v for k,v in inverted_homodimer_domains2.items() if k in unis }
    ALL_domains = readDomains('CATH-B' if periodic else 'CATH')
    chain_map = chainMap()

    count=0
    total=0
    for _, row in df.iterrows():
        pdb = row['PDB ID'] if periodic else row['id'][:4]
        domains = ALL_domains[pdb]

        ##empty domain, no point continuing
        

        ##get all subunits in heteromeric interactions and are in reduced topology

        interactions = {edge for pair, active in zip(row[key_interfaces].split(','),row[key_reduced].split(',')) if (pair[0]!=pair[2] and active == '1') for edge in pair.split('-') } if periodic else row['chains']


        for edge in interactions:
            total+=1
            ##relabel domains if needed
            if pdb in chain_map and edge in chain_map[pdb]:
                edge = chain_map[pdb][edge]
            if not domains:
                continue
            ##no information on this subunit for domains
            if edge not in domains:
                continue
            
            if tuple(domains[edge]) in inverted_homodimer_domains:
                count+=1
    return count,total

def visualiseOverlap(data_in):
    data = json.load(open(data_in,'r'))

    f, (ax1,ax2) = plt.subplots(2,1)
    #ax1t = ax1.twinx()

    total = list(data.values())[0][1]
    ax1.scatter(data.keys(), [v[0] for v in data.values()])
    #ax1t.scatter(data.keys(), [v[0]/total for v in data.values()])

    unis=getUniqueHomodimerDomains()
    s2 = invertDomains(readDomains('HCath'))
    #s = [list(getUniqueHomodimerDomains(i+1)-getUniqueHomodimerDomains(i))[0] for i in range(len(unis))]

    c = Counter(np.ediff1d([v[0] for v in data.values()]))

    #matched = {k:v for (k,v) in zip(s,np.ediff1d([v[0] for v in data.values()]))}
    #return matched
    ax2.scatter(sorted(c.keys()),[c[k] for k in sorted(c.keys())])
    ax2.set_yscale('log')
    plt.show(block=False)
    

def linkedProteinGenerator(df):
    key_interfaces = 'List of interface types (all identical subunits are given the same code)'
    key_reduced = 'Interfaces included in the reduced topology (1 - yes; 0 - no)'

    inverted_homodimer_domains = invertDomains(readDomains('HCath'))
    ALL_domains = readDomains('CATH-B')
    chain_map = chainMap()


    for _, row in df.iterrows():
        pdb = row['PDB ID']
        domains = ALL_domains[pdb]

        ##empty domain, no point continuing
        if not domains:
            continue

        ##get all subunits in heteromeric interactions and are in reduced topology
        interactions={edge for pair, active in zip(row[key_interfaces].split(','),row[key_reduced].split(',')) if (pair[0]!=pair[2] and active == '1') for edge in pair.split('-') }

        for edge in interactions:
            ##relabel domains if needed
            if pdb in chain_map and edge in chain_map[pdb]:
                edge = chain_map[pdb][edge]
                
            ##no information on this subunit for domains
            if edge not in domains:
                continue

            if tuple(domains[edge]) in inverted_homodimer_domains:
                for comp in inverted_homodimer_domains[tuple(domains[edge])]:
                    yield (pdb + '_' + edge, '{}_{}'.format(*comp.split('_')))
                        
    return

def randomProteinSampler(df, N_SAMPLE_LIMIT=100000, REQ_domain_info=True, ENFORCE_nonmatch=False):

    key_interfaces = 'List of interface types (all identical subunits are given the same code)'
    key_reduced = 'Interfaces included in the reduced topology (1 - yes; 0 - no)'

    inverted_homodimer_domains = invertDomains(readDomains('HCath'))
    homodimer_set = [pdb for proteins in inverted_homodimer_domains.values() for pdb in proteins]
    heterodimer_set = []
    heteromer_domains = readDomains('CATH-B')
    homodimer_domains = readDomains('HCath')
    chain_map = chainMap()

    
    for _, row in df.iterrows():
        pdb = row['PDB ID']
        domains = heteromer_domains[pdb]

        ##empty domain, no point continuing
        if REQ_domain_info and not domains:
            continue

        interactions={edge for pair, active in zip(row[key_interfaces].split(','),row[key_reduced].split(',')) if (pair[0]!=pair[2] and active == '1') for edge in pair.split('-')}
        
        for edge in interactions:
            ##relabel domain
            if pdb in chain_map and edge in chain_map[pdb]:
                edge = chain_map[pdb][edge]

                ##no information on this subunit for domains
                if REQ_domain_info and edge not in domains:
                    continue

                if not REQ_domain_info or tuple(domains[edge]) in inverted_homodimer_domains:
                    heterodimer_set.append('{}_{}'.format(pdb,edge))

    if ENFORCE_nonmatch:
        yield_count = 0
        while yield_count <= N_SAMPLE_LIMIT:
            het_option = choice(heterodimer_set)
            hom_option = choice(homodimer_set)
            while het_option[-1] in heteromer_domains[het_option[:4]] and hom_option[-1] in homodimer_domains[hom_option[:4]] and heteromer_domains[het_option[:4]][het_option[-1]] == homodimer_domains[hom_option[:4]][hom_option[-1]]:
                het_option = choice(heterodimer_set)
                hom_option = choice(homodimer_set)

            yield_count += 1
            yield (het_option, hom_option)
            
    else:
        for hetero_index, homo_index_ in zip(randint(0,len(heterodimer_set),N_SAMPLE_LIMIT),randint(0,len(homodimer_set),N_SAMPLE_LIMIT)):
            yield (heterodimer_set[hetero_index], homodimer_set[homo_index_])
    
    return

   

def chainMap():
    chainmap = defaultdict(dict)
    with open('chain_map.txt','r') as file_:
        for line in file_:
            (pdb, file_chain, pdb_chain) = line.rstrip().split('\t')
            chainmap[pdb][file_chain] = pdb_chain
    return dict(chainmap)


def main():

    parser = argparse.ArgumentParser(description = 'Domain comparison suite')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-P", "--parallel", action="store_true",dest='exec_style')
    group.add_argument("-S", "--sequential", action="store_false",dest='exec_style')

    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("-D", "--domain", action="store_true",dest='exec_mode')
    group2.add_argument("-R", "--random", action="store_false",dest='exec_mode')

    parser.add_argument('--R_mode', type=int,dest='R_mode')
    parser.add_argument('-N','--N_samples', type=int,dest='N_samples')
    parser.add_argument('--file_name', type=str,dest='file_name')
    parser.set_defaults(exec_style=False,exec_mode=True,N_samples=None,file_name=None,R_mode=0)
    args = parser.parse_args()


    df =readHeteromers()

    random_condition={0: (False,False), 1: (True,False), 2: (True,True), 3: (False,True)}
    proteinGenerator = linkedProteinGenerator(df) if args.exec_mode else randomProteinSampler(df,args.N_samples,*random_condition[args.R_mode])
    
    
    if args.exec_style:
        results = paralleliseAlignment(proteinGenerator)
    else:
        print('Running sequentially')
        results = {}
        for pdb in proteinGenerator:
            try:
                results = calculatePvalue(pdb)
            except Exception as err:
                print('Error on {}'.format(pdb),err)

            results['{}_{}_{}_{}'.format(*results[0])] = results[1]
  
    with open('{}_comparison.dict'.format(args.file_name or ('domain_match' if args.exec_mode else 'random')),'w') as f_out:
        json.dump(results,f_out)
        
if __name__ == "__main__":
    main()


