##local imports
from domains import readDomains, invertDomains, getUniqueHomodimerDomains,domainMatch
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
from itertools import combinations,product

def overlaps(df):
    rows = []

    for _,row in df.iterrows():
        for pair in row['interfaces']: #combinations(row['domains'].keys()):
            if isinstance(row['BSAs'],str):
                try:
                    row['BSAs'] = eval(row['BSAs'])
                except:
                    print(type(row['BSAs']))
                    print("ERR",row['BSAs'])
            #),['C-D'],pair,type(pair))
            #row['BSAs'][pair])#domainMatch(row['domains'],*pair.split('-')))
            try:
                rows.append({'id':row['PDB_id']+'_'+str(pair),'shared':domainMatch(row['domains'],*pair.split('-')),'BSA':row['BSAs'][pair]})
            except:
                print(row)

    print('concant')
    return pandas.DataFrame(rows)
            

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


def checkINTExisting(df_HET, df_HOM):
    for clean_df in (df_HET, df_HOM):
        drop_indices = [ind for ind,row in clean_df.iterrows() if not os.path.exists(f'/scratch/asl47/PDB/INT/{row["PDB_id"].upper()}.int')]
        print('Had to drop ',len(drop_indices))
        clean_df.drop(drop_indices,inplace=True)


def getDomainPWeights(domains,samples=None):
    lengths = np.array([len(v) for v in domains.values()])
    return lengths/np.sum(lengths)
    
    
def RPS_wrapper(df_HET, df_HOM, N_SAMPLE_LIMIT,match_partials=False):

    unique_comparisons = set()
    RPS = newRPS(df_HET, df_HOM, N_SAMPLE_LIMIT,match_partials)

    while len(unique_comparisons) < N_SAMPLE_LIMIT:
        unique_comparisons.add(next(RPS))
    return unique_comparisons
        

def newRPS(df_HET, df_HOM, N_SAMPLE_LIMIT,match_partials=False):
    inverted_domains_HET = invertCSVDomains(df_HET,match_partials)
    inverted_domains_HOM = invertCSVDomains(df_HOM,match_partials)

    weights_HET = getDomainPWeights(inverted_domains_HET)
    #weights_HOM = getDomainPWeights(inverted_domains_HOM)

    while True:
        sampled_domain = choice(list(inverted_domains_HET),p=weights_HET)
        if sampled_domain not in inverted_domains_HOM:
            continue
        
        chain_HET = choice(inverted_domains_HET[sampled_domain])
        chain_HOM = choice(inverted_domains_HOM[sampled_domain])
        
        if chain_HET == chain_HOM:
            #print('Self-compare, skip')
            continue
        yield (chain_HET,chain_HOM)
        
def newFPS(df_HET, df_HOM,match_partials=False):
    inverted_domains_HET = invertCSVDomains(df_HET,match_partials)
    inverted_domains_HOM = invertCSVDomains(df_HOM,match_partials)

    for sampled_domain in inverted_domains_HET:
        if sampled_domain not in inverted_domains_HOM:
            continue
        
        for (chain_HET,chain_HOM) in product(inverted_domains_HET[sampled_domain],inverted_domains_HOM[sampled_domain]):
        
            if chain_HET == chain_HOM:
                continue
            yield (chain_HET,chain_HOM)

def newEPS(df_HET, df_HOM, N_SAMPLE_LIMIT,match_partials=True):
    inverted_domains_HET = invertCSVDomains(df_HET,match_partials)
    inverted_domains_HOM = invertCSVDomains(df_HOM,match_partials)

    anti_HET = generateAntiDomains(df_HET)
    anti_HOM = generateAntiDomains(df_HOM)

    weights_HET = getDomainPWeights(inverted_domains_HET)
    weights_HOM = getDomainPWeights(inverted_domains_HOM)

    yielded_samples = 0
    while yielded_samples < N_SAMPLE_LIMIT:
        sampled_domain_HET = choice(list(inverted_domains_HET),p=weights_HET)
        chain_HET = choice(inverted_domains_HET[sampled_domain_HET])
        
        sampled_domain_HOM = choice(list(inverted_domains_HOM),p=weights_HOM)
        chain_HOM = choice(inverted_domains_HOM[sampled_domain_HOM])
        
        if chain_HET == chain_HOM:
            #print('Self-compare, skip')
            continue

        if not match_partials and sampled_domain_HET == sampled_domain_HOM:
            continue
        
        if match_partials and any(darch in sampled_domain_HOM for darch in anti_HET[chain_HET[:4]][chain_HET[5]]):
            continue
        
        yielded_samples += 1
        yield (chain_HET,chain_HOM)
        
    return

    
def generateAntiDomains(df):
    anti_domains = {}

    for _, row in df.iterrows():
        pdb = row['PDB_id']
        domains = row['domains']

        interactions = row['interfaces']
       
        antis = defaultdict(set)
        for interaction in interactions:
            for edge in interaction.split('-'):
                antis[edge] |= {arch for chain in interaction.split('-') for arch in domains[chain]}

        ##if the subunit has no domains, don't include ones it interacts with either
        anti_domains[pdb] = {k:v for k,v in antis.items() if v and k in domains}
    return anti_domains
        
    
    
def randomProteinSampler(df_HET, df_HOM, domain_mode, N_SAMPLE_LIMIT,match_partials=False):
    
    inverted_homodimer_domains = invertCSVDomains(df_HOM)
    
    yielded_samples = 0
    heteromer_domains, heteromer_anti_domains = set(), {}
    homodimer_set, heteromer_set = set(), set()

    def yieldOp(pdb,chain,comp):
        
        nonlocal yielded_samples
        print(yielded_samples)
        if domain_mode == 'match':
            print('yield')
            yielded_samples += 1
            yield (f'{pdb}_{chain}', comp)
            if N_SAMPLE_LIMIT and yielded_samples >= N_SAMPLE_LIMIT:
                print('RETURNING\n\nRETURNED')
                return
        else:
            heteromer_set.add(f'{pdb}_{chain}')
            homodimer_set.add(comp)

    for domains in df_HET['domains']:
        for archs in domains.values():
            if match_partials:
                heteromer_domains |= set(archs)
            heteromer_domains.add(archs)
            
    
    
    for _, row in df_HET.iterrows():
        pdb = row['PDB_id']
        domains = row['domains']

        ##empty domain, no point continuing
        #if not domains:
        #    continue

        interactions = row['interfaces']
        flat_chains = {chain for interaction in interactions for chain in interaction.split('-')}

        anti_domains = defaultdict(set)
        for interaction in interactions:
            for edge in interaction.split('-'):
                anti_domains[edge] |= {arch for chain in interaction.split('-') if chain in domains for arch in domains[chain]}

        ##if the subunit has no domains, don't include ones it interacts with either
        heteromer_anti_domains[pdb] = {k:v for k,v in anti_domains.items() if v and k in domains}
        
        
        #if domain_mode == 'match':
        for chain in flat_chains:
            ##no information on this subunit for domains
            #if chain not in domains:
            #    continue
                    
            if match_partials:
                for partial_domain in set(domains[chain]):
                    if (partial_domain,) in inverted_homodimer_domains:
                        for comp in inverted_homodimer_domains[(partial_domain,)]:
                            yield from yieldOp(pdb,chain,comp)
                                
            else:
                if domains[chain] in inverted_homodimer_domains:
                    for comp in inverted_homodimer_domains[domains[chain]]:
                        yield from yieldOp(pdb,chain,comp)
               
    if domain_mode == 'match':
        return
    
    #homodimer_set = list(f'{p}_{list(c)[0]}' for p,c,D in zip(homodimer_table['PDB_id'],homodimer_table['interfaces'],homodimer_table['domains']) if (list(c)[0] in D and D[list(c)[0]] in heteromer_domains))
    #heteromer_set = list(f'{p}_{j}' for p,c,D in zip(df['PDB_id'],df['interfaces'],df['domains']) for i in c for j in i.split('-') if j in D)

    homodimer_domains = {}
    for _,row in df_HOM.iterrows():
        if row['domains'] is not None:
            homodimer_domains[row['PDB_id']] = row['domains']

    homodimer_set, heteromer_set = (tuple(set_) for set_ in (homodimer_set, heteromer_set))
    
    if domain_mode == 'enforce':
        while yielded_samples < N_SAMPLE_LIMIT:
            het_option = choice(heteromer_set)
            hom_option = choice(homodimer_set)
            
            try:
                HADs = heteromer_anti_domains[het_option[:4]][het_option[-1]]
                BADs = homodimer_domains[hom_option[:4]][hom_option[-1]]
                if any(darch in HADs for darch in BADs):
                    continue
            except KeyError:
                continue # pass

            yield (het_option, hom_option)
            yielded_samples += 1
        return

    ##elif domain_mode =='random'
    else: 
        for pairing in zip(choice(heteromer_set,N_SAMPLE_LIMIT,True),choice(homodimer_set,N_SAMPLE_LIMIT,True)):
            yield pairing
        return  

def chainMap():
    chainmap = defaultdict(dict)
    with open('chain_map.txt','r') as file_:
        for line in file_:
            (pdb, file_chain, pdb_chain) = line.rstrip().split('\t')
            chainmap[pdb][file_chain] = pdb_chain
    return dict(chainmap)


def main(args):
    #if args.exec_source:
    #    print('Loading periodic data')
    #    df = loadCSV('Periodic_heteromers_C.csv')
    #else:
    #    print('Loading PDB heterodimer data')
    #    df = loadCSV('PDB_heterodimers.csv')

    df = loadCSV('New_heteromersR.csv')
    df2 = loadCSV('New_homomersR.csv')

    #checkINTExisting(df,df2)

    if args.exec_mode == 'match':
        if args.N_samples is None:
            print('Exhaustive iterating domain matched comparisons')
            proteinGenerator = newFPS(df,df2,args.allow_partials)
        else:
            print('Sampling domain matched comparisons')
            proteinGenerator = RPS_wrapper(df,df2,args.N_samples,args.allow_partials)
    else:
        print('Sampling anti-domain enforced comparisons')
        proteinGenerator = newEPS(df,df2,args.N_samples,args.allow_partials)
        
    if args.exec_style:
        results = paralleliseAlignment(proteinGenerator)
    else:
        print('Running sequentially')
        results = {}
        for pdb in proteinGenerator:
            try:
                single_result = calculatePvalue(pdb)
                os.remove('/scratch/asl47/PDB/NEEDLE/{}_{}_{}_{}.needle'.format(*single_result[0]))
            except Exception as err:
                print('Error on {}'.format(pdb),err)
            

            results['{}_{}_{}_{}'.format(*single_result[0])] = single_result[1]
  

    if args.json:
        with open('{}_{}_comparison.dict'.format('Table' if args.exec_source else 'PDB', args.file_name or ('domain_match' if args.exec_mode else 'random')),'w') as f_out:
            json.dump(results,f_out)
    else:
        with open(f'{args.file_name}_comparison.csv','w') as f_out:
            df = pandas.DataFrame(results)
            df.columns = ['pval_F','pval_S','pval_T','hits','similarity','score','align_length','overlap']
            df.to_csv(f_out,index=False,columns=['pval_F','pval_S','pval_T','hits','similarity','score','align_length','overlap'])
            #np.array([val for val in results.values() if val!='error'],dtype=float).tofile(f_out)
        
        
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
    parser.add_argument('--json', type=bool,dest='json')
    parser.add_argument('-N','--N_samples', type=int,dest='N_samples')
    parser.add_argument('--file_name', type=str,dest='file_name')
    parser.add_argument('--partial', action='store_true',dest='allow_partials')
    parser.set_defaults(exec_style=False,exec_mode=None,exec_source=True,N_samples=None,file_name=None,allow_partials=False,json=False)
    
    args = parser.parse_args()

    if not args.exec_mode:
        print('Defaulting to random alignment')
        args.exec_mode = 'random'
 
    if args.exec_mode != 'match' and not args.N_samples:
        print('Random sampling amount defaulted to 10,000')
        args.N_samples = 10000


    main(args)


