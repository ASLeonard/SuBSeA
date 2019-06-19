import pandas
from domains import readDomains, invertDomains
from SubSeA import paralleliseAlignment, calculatePvalue
from numpy.random import randint
from collections import defaultdict
import sys
import json
import argparse

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
    
    #parser.add_argument('-','--N_samples',default=None, type=int,help='take N random samples')
    parser.add_argument('-N','--N_samples',default=None, type=int)
    parser.set_defaults(exec_style=False,exec_mode=True)
    args = parser.parse_args()


    df =readHeteromers()

    proteinGenerator = linkedProteinGenerator(df) if args.exec_mode else randomProteinSampler(df,args.N_samples)
    
    
    if args.exec_style:
        dic = paralleliseAlignment(proteinGenerator)
    else:
        print('Running sequentially')
        dic = {}
        for pdb in proteinGenerator:
            try:
                results = calculatePvalue(pdb)
            except Exception as err:
                print('Error on {}'.format(pdb),err)

            dic['{}_{}_{}_{}'.format(*results[0])] = results[1]
  
    with open('{}_run_1.dict'.format('domain' if args.exec_mode else 'random'),'w') as f_out:
        json.dump(dic,f_out)
        
if __name__ == "__main__":
    main()


