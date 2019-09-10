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

def scrapePDBs(df):
    with open('period_pdb_codes.txt', 'w') as file_out:
          file_out.write(', '.join(df['PDB ID']))



def getDomainPWeights(domains,samples=None):
    lengths = np.array([len(v) for v in domains.values()])
    return lengths/np.sum(lengths)
    
    
def RPS_wrapper(df_HET, df_HOM, N_SAMPLE_LIMIT,match_partials=False):

    unique_comparisons = set()
    RPS = newRPS(df_HET, df_HOM, match_partials)

    while len(unique_comparisons) < N_SAMPLE_LIMIT:
        unique_comparisons.add(next(RPS))
    return unique_comparisons
        

def newRPS(df_HET, df_HOM, match_partials=False):
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
        
        if chain_HET[:4] == chain_HOM[:4]:
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
        
            if chain_HET[:4] == chain_HOM[:4]:
                continue
            yield (chain_HET,chain_HOM)

##dead
from scipy.stats import chi2_contingency
def domainSubunitRewire(df):

    interaction_edges = []

    for _,row in df.iterrows():
        domains = row['domains']
        interactions = row['interfaces']
        for interaction in interactions:
            interaction_edges.extend([domains[C] for C in interaction.split('-')])
    print(len(interaction_edges))
    
    
    #chi2_contingency(obs, lambda_="log-likelihood")
    return 'dead'

from collections import Counter
def duplicateIntersection(domain_1, domain_2):
    overlap = Counter(domain_1) & Counter(domain_2)
    return tuple(dom for dom in sorted(overlap.elements()))

def conv(dq):
    a=[]
    for _,row in dq.iterrows():
        for unit in row['interfaces'].split('-'):
            if unit not in row['domains']:
                continue
        #try:
        #    domains = {K[0]:K[2:] for K in row['domains'].split(';')}
        #except:
        #    continue
    
        a.append({'PDB_id':row['PDB_id'],'interfaces':[row['interfaces']],'domains':row['domains'],'BSAs':{row['interfaces']:row['BSAs']}})

    return pandas.DataFrame(a)

def sharedMatchingAlgorithm2(df_HET, df_HOM):
    inverted_domains_HOM_full = invertCSVDomains(df_HOM,False,True)
    inverted_domains_HOM_partial = invertCSVDomains(df_HOM,True,True)

    
    comparisons_to_make = []
    fail_c=0
    for _, row in df_HET.iterrows():

        
        pdb = row['PDB_id']
        domains = row['domains']
        interactions = row['interfaces']

        for interaction_pair in interactions:
            subunits = interaction_pair.split('-')
            mutual_domains = duplicateIntersection(*(domains[C] for C in subunits))
            for S1, S2 in zip(subunits,reversed(subunits)):
                if mutual_domains == domains[S1]:
                    
                    if mutual_domains in inverted_domains_HOM_full:
                        for comp_pdb in inverted_domains_HOM_full[mutual_domains]:
                            if pdb == comp_pdb[:4]:
                                continue
                            
                            comparisons_to_make.append((f'{pdb}_{S1}_{S2}',comp_pdb,'MUT'))
                    else:
                        fail_c+=1
                elif mutual_domains:
                    if mutual_domains in inverted_domains_HOM_full:
                        for comp_pdb in inverted_domains_HOM_full[mutual_domains]:
                            if pdb == comp_pdb[:4]:
                                continue
                            
                            comparisons_to_make.append((f'{pdb}_{S1}_{S2}',comp_pdb,'MPA'))
                    
                else:
                     if domains[S2] in inverted_domains_HOM_full:
                        for comp_pdb in inverted_domains_HOM_full[domains[S2]]:
                            if pdb == comp_pdb[:4]:
                                continue
                            comparisons_to_make.append((f'{pdb}_{S1}_{S2}',comp_pdb,'DNO'))

    print(fail_c)
    return comparisons_to_make
                    

def filterDataset(df,thresh,hom_mode=False):
    if thresh not in {50,70,90}:
        print('not implemented')
        return
    
    with open(f'PDB_clusters_{thresh}.txt') as cluster_file:
        redundant_pdbs = [set(line.split()) for line in cluster_file]


    used_cluster_interactions = set()
    new_df = []
    
    for _,row in df.iterrows():
        pdb = row['PDB_id']
        unique_interactions = []

        
        for interaction_pair in row['interfaces']:
            cluster_indexes = []
            for chain in interaction_pair.split('-'):
                for index, cluster in enumerate(redundant_pdbs):
                    if f'{pdb.upper()}_{chain}' in cluster:
                        cluster_indexes.append(index)
                        break

            cluster_indexes = tuple(cluster_indexes)
            if hom_mode and cluster_indexes[0]!=cluster_indexes[1]:
                print('not right',pdb, interaction_pair)
            
            if cluster_indexes not in used_cluster_interactions:
                used_cluster_interactions.add(cluster_indexes)
                unique_interactions.append(interaction_pair)

        if unique_interactions:
            flat_chains = {C for pair in unique_interactions for C in pair.split('-')}
            new_df.append({'PDB_id':pdb,'interfaces':unique_interactions,'BSAs':{pair:row['BSAs'][pair] for pair in unique_interactions},'domains':{C:row['domains'][C] for C in flat_chains}})

    return pandas.DataFrame(new_df)
    
    

    

def sharedMatchingAlgorithm(df_HET, df_HOM):
    inverted_domains_HOM_full = invertCSVDomains(df_HOM)
    inverted_domains_HOM_partial = invertCSVDomains(df_HOM,True)
    
    comparisons_to_make = []
    ss=0
    for _, row in df_HET.iterrows():

        mutual_comparisons, matched_comparisons = set(), set()
        partial_comparisons_S, full_comparisons_S = set(), set()
        partial_comparisons_N, full_comparisons_N = set(), set()

        
        pdb = row['PDB_id']
        domains = row['domains']

        interactions = row['interfaces']
        flat_chains = {chain for interaction in interactions for chain in interaction.split('-')}

        partner_domains = defaultdict(dict)

        shared=[]
        for flat_chain in flat_chains:
            chain_domains = domains[flat_chain]
            ## mutual overlaps
            
            for interaction in interactions:
                if flat_chain not in interaction:
                    continue
                for edge in interaction.split('-'):
                    if edge == flat_chain:
                        continue

                    mutual_domains = duplicateIntersection(chain_domains,domains[edge])
                    if mutual_domains:
                        shared.append((flat_chain,edge))
                        if mutual_domains in inverted_domains_HOM_full:
                            for pdb_hom in inverted_domains_HOM_full[mutual_domains]:
                                if pdb_hom[:4] != pdb:
                                    if len(mutual_domains) == len(chain_domains):
                                        matched_comparisons.add((f'{pdb}_{flat_chain}',pdb_hom))
                                    else:
                                        mutual_comparisons.add((f'{pdb}_{flat_chain}',pdb_hom))
                    else:
                        if domains[edge] in inverted_domains_HOM_full:
                            for pdb_hom in inverted_domains_HOM_full[domains[edge]]:
                                if pdb_hom[:4] != pdb:
                                    full_comparisons_N.add((f'{pdb}_{flat_chain}',pdb_hom))
                        if len(domains[edge]) >= 2:
                            for domain in domains[edge]:
                                if domain in inverted_domains_HOM_partial:
                                    for pdb_hom in inverted_domains_HOM_partial[domain]:
                                        if pdb_hom[:4] != pdb and (f'{pdb}_{flat_chain}',pdb_hom) not in full_comparisons_N:
                                            partial_comparisons_N.add((f'{pdb}_{flat_chain}',pdb_hom))
                                                           

            ## per full chain domain comparisons

            if chain_domains in inverted_domains_HOM_full:
                for pdb_hom in inverted_domains_HOM_full[chain_domains]:
                    if pdb_hom[:4] == pdb:
                        continue
                    comp_tuple = (f'{pdb}_{flat_chain}',pdb_hom)
                    if comp_tuple not in matched_comparisons:
                        full_comparisons_S.add(comp_tuple)

            ## partial comparisons
            if len(chain_domains) >= 2:
                for chain_domain in chain_domains:
                    if chain_domain in inverted_domains_HOM_partial:
                        for pdb_hom in inverted_domains_HOM_partial[chain_domain]:
                            if pdb_hom[:4] == pdb:
                                continue
                            comp_tuple = (f'{pdb}_{flat_chain}',pdb_hom)
                            if comp_tuple not in full_comparisons_S and comp_tuple not in mutual_comparisons and comp_tuple not in matched_comparisons:
                                partial_comparisons_S.add(comp_tuple)
                    
        #print(pdb,shared)
        if shared:
            ss+=1
        #filter out mistakes
        full_comparisons_N -= full_comparisons_S
        partial_comparisons_N -= partial_comparisons_S
        partial_comparisons_N -= full_comparisons_N
        
        for data, code in zip((matched_comparisons,mutual_comparisons,full_comparisons_S,full_comparisons_N,partial_comparisons_S,partial_comparisons_N),('MF','MP','FS','FN','PS','PN')):
            comparisons_to_make.extend([comp+(code,) for comp in data])
    print('Shared domains on ',ss)
    return comparisons_to_make
        
def doubleEntries(comps,sorted_mode=False):
    codes= set()
    types= {}
    for comp in comps:
        val = tuple(sorted((comp[0],comp[1]))) if sorted_mode else (comp[0],comp[1])
        if val not in codes:
            codes.add(val)
            types[val]=comp[2]
        elif types[val]!=comp[2]:
            print(types[val],comp[2],val)

def getUniqueChains(df):
    pdb_chains = []

    for _,row in df.iterrows():
        interactions = row['interfaces']
        chains = {chain for interaction in interactions for chain in interaction.split('-')}
        pdb_chains.extend([f"{row['PDB_id']}_{chain}" for chain in chains])

    return pdb_chains
    
def EPS_wrapper(df_HET, df_HOM, N_SAMPLE_LIMIT,):

    unique_comparisons = set()
    EPS = newEPS(df_HET, df_HOM)

    while len(unique_comparisons) < N_SAMPLE_LIMIT:
        unique_comparisons.add(next(EPS))
    return unique_comparisons

def newEPS(df_HET, df_HOM):
    het_options, hom_options = getUniqueChains(df_HET), getUniqueChains(df_HOM)

    anti_domains = generateAntiDomains(df_HET,df_HOM)
    while True:

        chain_HET = choice(het_options)
        chain_HOM = choice(hom_options)
        
        if chain_HET[:4] == chain_HOM[:4]:
            #print('Self-compare, skip')
            continue

        if any(arch in anti_domains[chain_HOM[:4]][chain_HOM[5]] for arch in anti_domains[chain_HET[:4]][chain_HET[5]]):
            continue

        yield (chain_HET,chain_HOM)

    
def generateAntiDomains(df,df2):
    anti_domains = {}

    for pdb in set(df['PDB_id']).union(set(df2['PDB_id'])):
        domains, interactions = {}, []

        for frame in (df,df2):
            data = frame.loc[frame['PDB_id']==pdb]
            if not data.empty:
                domains = {**domains, **data.iloc[0]['domains']}
                interactions.extend(data.iloc[0]['interfaces'])
       
        antis = defaultdict(set)
        for interaction in interactions:
            for edge in interaction.split('-'):
                antis[edge] |= {arch for chain in interaction.split('-') for arch in domains[chain]}

        ##if the subunit has no domains, don't include ones it interacts with either
        anti_domains[pdb] = dict(antis)#{k:v for k,v in antis.items() if v and k in domains}
    return anti_domains

def shuffledInteractionDomains(df):

    table_of_observations = [[0,0],[0,0]]
    interaction_edges = []

    for domains, interactions in zip(df['domains'],df['interfaces']):
        for interaction in interactions:
            local_domains = [domains[C] for C in interaction.split('-')]
            table_of_observations[0][duplicateIntersection(*local_domains) != ()] += 1
            interaction_edges.extend(local_domains)

    interaction_edges = np.asarray(interaction_edges)
    ##shuffle it
    np.random.shuffle(interaction_edges)
    interaction_edges = interaction_edges.reshape((2,-1))

    overlap_results = np.apply_along_axis(lambda x: duplicateIntersection(*x)!=(),arr=interaction_edges,axis=0)
    for form in (0,1):
        table_of_observations[1][form] = len(overlap_results[overlap_results==form])
    
    

    return table_of_observations, fisher_exact(table_of_observations)

from scipy.stats import fisher_exact
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

    assert args.filter_level in {50,70,90,100}, 'Invalid filter level'
    if args.filter_level == 100:
        args.filter_level = 'unfiltered'
    
    df = loadCSV(f'Heteromers_{args.filter_level}.csv')
    df2 = loadCSV(f'Homomers_{args.filter_level}.csv')


    if args.exec_mode == 'match':
        proteinGenerator = sharedMatchingAlgorithm2(df,df2)
        if args.N_samples:
            proteinGenerator = [proteinGenerator[index] for index in choice(len(proteinGenerator),replace=False,size=args.N_samples)]
        #if args.N_samples is None:
        #    print('Exhaustive iterating domain matched comparisons')
        #    proteinGenerator = newFPS(df,df2,args.allow_partials)
        #else:
        #    print('Sampling domain matched comparisons')
        #    proteinGenerator = RPS_wrapper(df,df2,args.N_samples,args.allow_partials)
    else:
        print('Sampling anti-domain enforced comparisons')
        proteinGenerator = EPS_wrapper(df,df2,args.N_samples)
        
    if args.exec_style:
        results = paralleliseAlignment(proteinGenerator,args.file_name)
    else:
        print('Running sequentially')
        results = []
        for pdb in proteinGenerator:
            try:
                single_result = calculatePvalue(pdb)
                os.remove('/scratch/asl47/PDB/NEEDLE/{}_{}_{}_{}.needle'.format(*single_result[0]))
            except Exception as err:
                print('Error on {}'.format(pdb),err)
            

            results.append((single_result[0],)+single_result[1])
  

    if args.json:
        with open('{}_{}_comparison.dict'.format('Table' if args.exec_source else 'PDB', args.file_name or ('domain_match' if args.exec_mode else 'random')),'w') as f_out:
            json.dump(results,f_out)
    else:

        columns = ['id','code','pval_F','pval_S','pval_T','pval_F2','pval_S2','pval_T2','hits','similarity','score','align_length','overlap']
        df = pandas.DataFrame(results)
        df.columns = columns
        
        with open(f'/rscratch/asl47/PDB_results/{args.file_name}_comparison.csv','w') as f_out:
            
            df.to_csv(f_out,index=False,columns=columns)

        
        
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
    parser.add_argument('--filter', type=int,dest='filter_level')
    parser.add_argument('-N','--N_samples', type=int,dest='N_samples')
    parser.add_argument('--file_name', type=str,dest='file_name')
    parser.add_argument('--partial', action='store_true',dest='allow_partials')
    parser.set_defaults(exec_style=False,exec_mode=None,exec_source=True,N_samples=None,file_name=None,allow_partials=False,json=False,filter_level=50)
    
    args = parser.parse_args()

    if not args.exec_mode:
        print('Defaulting to random alignment')
        args.exec_mode = 'random'
 
    if args.exec_mode != 'match' and not args.N_samples:
        print('Random sampling amount defaulted to 10,000')
        args.N_samples = 10000


    main(args)


