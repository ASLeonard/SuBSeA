##local imports
from domains import readDomains, invertDomains, getUniqueHomodimerDomains
from SubSeA import paralleliseAlignment, calculatePvalue
from pdb_visualisation import scrubInput

##global imports
import pandas
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

def randomProteinSampler(df, periodic_table=True, N_SAMPLE_LIMIT=100000, REQ_domain_info=True, ENFORCE_nonmatch=False):

    key_interfaces = 'List of interface types (all identical subunits are given the same code)'
    key_reduced = 'Interfaces included in the reduced topology (1 - yes; 0 - no)'

    
    heteromer_domains = readDomains('CATH-B' if periodic_table else 'CATH')
    homodimer_domains = readDomains('HCath')
    inverted_homodimer_domains = invertDomains(homodimer_domains)
    homodimer_set = [pdb for proteins in inverted_homodimer_domains.values() for pdb in proteins]# if REQ_domain_info else list(homodimer_domains.keys())
    heterodimer_set = []
    chain_map = chainMap()

    
    for _, row in df.iterrows():
        pdb = row['PDB ID'] if periodic_table else row['id'][:4]
        domains = heteromer_domains[pdb]

        ##empty domain, no point continuing
        if REQ_domain_info and not domains:
            continue

        interactions={edge for pair, active in zip(row[key_interfaces].split(','),row[key_reduced].split(',')) if (pair[0]!=pair[2] and active == '1') for edge in pair.split('-')} if periodic_table else row['chains']
        
        for edge in interactions:
            ##relabel domain
            if pdb in chain_map and edge in chain_map[pdb]:
                edge = chain_map[pdb][edge]

            ##no information on this subunit for domains
            if REQ_domain_info:
                if edge in domains and tuple(domains[edge]) in inverted_homodimer_domains:
                    if N_SAMPLE_LIMIT<=0:
                        for comp in inverted_homodimer_domains[tuple(domains[edge])]:
                            yield (pdb + '_' + edge, '{}_{}'.format(*comp.split('_')))
                    else:
                        heterodimer_set.append('{}_{}'.format(pdb,edge))
            else:
                heterodimer_set.append('{}_{}'.format(pdb,edge))

    if N_SAMPLE_LIMIT <= 0:
        return
    elif ENFORCE_nonmatch:
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


def main(args):
    if args.exec_source:
        print('Loading periodic data')
        df = readHeteromers()
    else:
        print('Loading PDB heterodimer data')
        df = scrubInput('PDB_Heterodimers',0,1,'CATH')
    
        
    random_condition = {0: (False,False), 1: (True,False), 2: (True,True), 3: (False,True)}
    
    proteinGenerator = randomProteinSampler(df,args.exec_source,args.N_samples,*random_condition[args.R_mode])
    
    
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
    group2.add_argument("-D", "--domain", action="store_true",dest='exec_mode')
    group2.add_argument("-R", "--random", action="store_false",dest='exec_mode')    

    parser.add_argument('--R_mode', type=int,dest='R_mode')
    parser.add_argument('-N','--N_samples', type=int,dest='N_samples')
    parser.add_argument('--file_name', type=str,dest='file_name')
    parser.set_defaults(exec_style=False,exec_mode=True,exec_source=True,N_samples=None,file_name=None,R_mode=-1)
    
    args = parser.parse_args()

    if args.exec_mode:
        args.N_samples = -1
        args.R_mode = 1
    elif not args.N_samples:
        print('Random sampling amount defaulted to 10,000')
        args.N_samples = 10000
        

    if args.R_mode == -1:
        args.R_mode = 0
        print('Random style unset, running as fully random')
        
    main(args)
