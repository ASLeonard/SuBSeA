##local imports
from domains import readDomains, invertDomains, getUniqueHomodimerDomains
from SubSeA import paralleliseAlignment, calculatePvalue
#from pdb_visualisation import scrubInput
from utility import loadCSV, invertCSVDomains

##global imports
import pandas
from numpy.random import randint, choice
from collections import defaultdict, Counter
import sys
import os
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

def heteromerHomodimerOverlap(df,periodic=True,H_max=None):
    key_interfaces = 'List of interface types (all identical subunits are given the same code)'
    key_reduced = 'Interfaces included in the reduced topology (1 - yes; 0 - no)'

    inverted_homodimer_domains = invertDomains(readDomains('HCath'))
    #unique_archs, unique_counts=getUniqueHomodimerDomains(10) 
    #inverted_homodimer_domains = {k:v for k,v in inverted_homodimer_domains2.items() if k in unis}
    ALL_domains = readDomains('CATH-B' if periodic else 'CATH')
    chain_map = chainMap()

    counts,uc = [], []
    #fs=[]
    
    for (unique_archs, unique_counts) in getUniqueHomodimerDomains(H_max):
        if unique_archs not in inverted_homodimer_domains:
            continue
        uc.append(unique_counts)
        count, total = 0,0

        for _, row in df.iterrows():
            pdb = row['PDB ID'] if periodic else row['id'][:4]
            domains = ALL_domains[pdb]

            interactions = {edge for pair, active in zip(row[key_interfaces].split(','),row[key_reduced].split(',')) if (pair[0]!=pair[2] and active == '1') for edge in pair.split('-') } if periodic else row['chains']


            for edge in interactions:
                total+=1
                ##return elabel domains if needed
                if pdb in chain_map and edge in chain_map[pdb]:
                    edge = chain_map[pdb][edge]
                if not domains:
                    continue
                ##no information on this subunit for domains
                if edge not in domains:
                    #fs.append('{}_{}'.format(pdb.upper(),edge))
                    continue
            
                if tuple(domains[edge]) == unique_archs:
                    count+=1
        #return fs
        counts.append(count)
        
    return {'matches':counts,'h_count':uc,'total':total}

def visualiseOverlap(data_in):
    data_raw = json.load(open(data_in,'r'))


    f, (ax1,ax2) = plt.subplots(2,1)
    

    total = data_raw['total']
    
    ax1.scatter(range(len(data_raw['matches'])),np.cumsum(data_raw['matches']),c=np.log10(data_raw['h_count']),cmap='viridis')
    ax1.set_xticks([])
    ax1.set_xlabel(r'$\leftarrow$ most common homodimer domain')
    ax1.set_ylabel('# matched heteromeric subunits')
    for cut in (.5,.75,.95):
        ax1.text(np.where(h>cut*np.sum(h))[0][0], 0, str(cut), fontsize=12)
    ax1t = ax1.twinx()
    ax1t.plot(np.cumsum(data_raw['matches'])/data_raw['total'],ls='--')
    ax1t.set_xticks([])
    #ax1t.set_ylims([0,np.cumsum(data_raw['matches'])])
    ax1t.set_ylabel('fraction of all matches')


    c = Counter(data_raw['matches'])
    cols=np.log10([sum(data_raw['h_count'][i] for i,j in enumerate(data_raw['matches']) if j==key) for key in sorted(c.keys())])
    
    ax2.scatter(sorted(c.keys()),[c[k] for k in sorted(c.keys())],c=cols,cmap='plasma',edgecolors='k')
    ax2.set_yscale('log')
    ax2.set_xscale('symlog')
    ax2.set_xlabel('# matched heteromeric subunits')
    ax2.set_ylabel('frequency')
    plt.show(block=False)
    

def randomProteinSampler(df, domain_mode, N_SAMPLE_LIMIT,match_partials=False):
    homodimer_table = loadCSV('PDB_homodimers.csv')

    for clean_df in (df, homodimer_table):
        drop_indices = [ind for ind,row in clean_df.iterrows() if not os.path.exists(f'/scratch/asl47/PDB/INT/{row["PDB_id"].upper()}.int')]
        clean_df.drop(drop_indices,inplace=True)
    
    inverted_homodimer_domains = invertCSVDomains(homodimer_table)
    
    yielded_samples = 0
    heteromer_domains, heteromer_anti_domains = {}, {}
    
    for _, row in df.iterrows():
        pdb = row['PDB_id']
        domains = row['domains']

        ##empty domain, no point continuing
        if not domains:
            continue

        interactions = row['interfaces']
        flat_chains = {chain for interaction in interactions for chain in interaction.split('-')}

        anti_domains = defaultdict(set)
        for chain in flat_chains:
            for interaction in interactions:
                if chain in interaction:
                    for edge in interaction.split('-'):
                        if edge in domains:
                            anti_domains[chain].add(domains[edge])

        heteromer_anti_domains[pdb] = dict(anti_domains)
        
        if domain_mode == 'match':
            for chain in flat_chains:
                ##no information on this subunit for domains
                if chain not in domains:
                    continue
                    
                if match_partials:
                    for partial_domain in domains[chain]:
                        if partial_domain in inverted_homodimer_domains:
                            for comp in inverted_homodimer_domains[partial_domain]:
                                yield (pdb + '_' + chain, '{}_{}'.format(*comp.split('_')))
                                yielded_samples += 1
                                if N_SAMPLE_LIMIT and yielded_samples > N_SAMPLE_LIMIT:
                                    return
                if chain in domains and domains[chain] in inverted_homodimer_domains:
                    for comp in inverted_homodimer_domains[domains[chain]]:
                        yield (pdb + '_' + chain, '{}_{}'.format(*comp.split('_')))
                        yielded_samples += 1
                        if N_SAMPLE_LIMIT and yielded_samples > N_SAMPLE_LIMIT:
                            return

    yield heteromer_anti_domains                
    if domain_mode == 'match':
        return
    
    homodimer_set = list(f'{p}_{list(c)[0]}' for p,c in zip(homodimer_table['PDB_id'],homodimer_table['interfaces']))
    heteromer_set = list(f'{p}_{j}' for p,c in zip(df['PDB_id'],df['interfaces']) for i in c for j in i.split('-'))

    homodimer_domains = {}
    for _,row in homodimer_table.iterrows():
        if row['domains'] is not None:
            homodimer_domains[row['PDB_id']] = row['domains']
    
    if domain_mode == 'enforce':
        while yielded_samples <= N_SAMPLE_LIMIT:
            het_option = choice(heteromer_set)
            hom_option = choice(homodimer_set)
            
            try:
                if any(darch in homodimer_domains[hom_option[:4]][hom_option[-1]] for darch in heteromer_anti_domains[het_option[:4]][het_option[-1]]):
                    continue
            except KeyError:
                pass

            yield (het_option, hom_option)
            yielded_samples += 1
        return
    
    else: #random domain_mode
        for hetero_index, homo_index_ in zip(randint(0,len(heteromer_set),N_SAMPLE_LIMIT),randint(0,len(homodimer_set),N_SAMPLE_LIMIT)):
            yield (heteromer_set[hetero_index], homodimer_set[homo_index_])
        return  

def chainMap():
    chainmap = defaultdict(dict)
    with open('chain_map.txt','r') as file_:
        for line in file_:
            (pdb, file_chain, pdb_chain) = line.rstrip().split('\t')
            chainmap[pdb][file_chain] = pdb_chain
    return dict(chainmap)


def main(args):
    if args.exec_source:
        print('Loading periodic data')
        df = loadCSV('Periodic_partial.csv')
        #loadCSV('Periodic_heteromers_CC.csv')#readHeteromers()
    else:
        print('Loading PDB heterodimer data')
        df = loadCSV('PDB_heterodimers.csv')

    
    proteinGenerator = randomProteinSampler(df,args.exec_mode,args.N_samples)
        
    if args.exec_style:
        results = paralleliseAlignment(proteinGenerator)
    else:
        print('Running sequentially')
        results = {}
        for pdb in proteinGenerator:
            #print(pdb)
            #try:
            single_result = calculatePvalue(pdb)
            #except Exception as err:
            #print('Error on {}'.format(pdb),err)

            results['{}_{}_{}_{}'.format(*single_result[0])] = single_result[1]
  
    with open('{}_{}_comparison.dict'.format('Table' if args.exec_source else 'PDB', args.file_name or ('domain_match' if args.exec_mode else 'random')),'w') as f_out:
        json.dump(results,f_out)

        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Domain comparison suite')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-P", "--parallel", action="store_true",dest='exec_style')
    group.add_argument("-S", "--sequential", action="store_false",dest='exec_style')

    group3 = parser.add_mutually_exclusive_group()
    group3.add_argument("-T", "--table", action="store_true",dest='exec_source')
    group3.add_argument("-H", "--heterodimers", action="store_false",dest='exec_source')
    
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("-D", "--domain", action="store_const",dest='exec_mode',const='match')
    group2.add_argument("-R", "--random", action="store_const",dest='exec_mode',const='random')
    group2.add_argument('-E','--enforce', action='store_const',dest='exec_mode',const='enforce')

    #parser.add_argument('--R_mode', type=int,dest='R_mode')
    parser.add_argument('-N','--N_samples', type=int,dest='N_samples')
    parser.add_argument('--file_name', type=str,dest='file_name')
    parser.set_defaults(exec_style=False,exec_mode=None,exec_source=True,N_samples=None,file_name=None,R_mode=-1)
    
    args = parser.parse_args()

    if not args.exec_mode:
        print('Defaulting to random alignment')
        args.exec_mode = 'random'
 
    if args.exec_mode != 'match' and not args.N_samples:
        print('Random sampling amount defaulted to 10,000')
        args.N_samples = 10000


    main(args)
