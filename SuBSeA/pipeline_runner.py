##local imports
from domains import duplicateIntersection
from binding_alignment import calculatePvalue
from utility import loadCSV, invertCSVDomains, makeDatasets, VALID_CLUSTERS
from pisa_XML import pullXML

##global imports
from functools import partial
from scipy.stats import fisher_exact
from multiprocessing import Pool
import os
import csv
import argparse
import numpy as np

def heteromericInteractionRunner(df_het):
    for _, row in df_het.iterrows():
        domains = row['domains']
        for interaction_pair in row['interfaces']:
            subunits = interaction_pair.split('-')
            mutual_domains = duplicateIntersection(*(domains[C] for C in subunits))
            code = ('MUT' if domains[subunits[0]]== domains[subunits[1]] else 'MPA') if mutual_domains else 'DNO'
            yield (f'{row["PDB_id"]}_{subunits[0]}_{subunits[1]}',f'{row["PDB_id"]}_{subunits[1]}_{subunits[0]}',code)
    return

def precursorInteractionRunner(df_HET, df_HOM):
    inverted_domains_HOM = invertCSVDomains(df_HOM,True,True)

    for _, row in df_HET.iterrows():
        pdb = row['PDB_id']
        domains = row['domains']
        interactions = row['interfaces']

        for interaction_pair in interactions:
            subunits = interaction_pair.split('-')
            mutual_domains = duplicateIntersection(*(domains[C] for C in subunits))
            code = ('MUT' if domains[subunits[0]]== domains[subunits[1]] else 'MPA') if mutual_domains else 'DNO'
            for subunit in subunits:
                if domains[subunit] in inverted_domains_HOM:
                    for comp_pdb in inverted_domains_HOM[domains[subunit]]:
                        if pdb != comp_pdb[:4]:
                            yield (f'{pdb}_{S1}_{S2}',comp_pdb,code)
    return

def shuffledInteractionDomains(df):
    table_of_observations = [[0,0],[0,0]]
    interaction_edges = []

    for domains, interactions in zip(df['domains'],df['interfaces']):
        for interaction in interactions:
            local_domains = [domains[C] for C in interaction.split('-')]
            table_of_observations[0][duplicateIntersection(*local_domains) != ()] += 1
            interaction_edges.extend(local_domains)
    interaction_edges = np.asarray(interaction_edges)

    np.random.shuffle(interaction_edges)
    interaction_edges = interaction_edges.reshape((2,-1))

    overlap_results = np.apply_along_axis(lambda x: duplicateIntersection(*x)!=(),arr=interaction_edges,axis=0)
    for form in (0,1):
        table_of_observations[1][form] = len(overlap_results[overlap_results==form])

    return table_of_observations, fisher_exact(table_of_observations)

def heteromericOverlapStats(df,df2):
    homomeric_domains = set()

    for domains in df2.domains:
        homomeric_domains |= {arch for archs in domains.values() for arch in archs}


    fractions = [[0,0],[0,0]]
    for _,row in df.iterrows():

        for interaction in row.interfaces:
            local_domains = [row.domains[C] for C in interaction.split('-')]
            ##HH type interaction
            overlap_domains = set(local_domains[0]) & set(local_domains[1])
            unique_domains = set(local_domains[0]) ^ set(local_domains[1])

            if overlap_domains:
                fractions[0][any(od in homomeric_domains for od in set(overlap_domains))]+=1
            else:
                fractions[1][any(od in homomeric_domains for od in set(unique_domains))]+=1
                
                
    print(f'Domain co-occurence\nHH: {fractions[0][0]}/{sum(fractions[0])} ({fractions[0][0]/sum(fractions[0]):.3f}%) \
            \nDH: {fractions[1][0]}/{sum(fractions[1])} ({fractions[1][0]/sum(fractions[1]):.3f}%)')
    return fractions

def paralleliseAlignment(pdb_pairs,file_name,domain,cluster):
    data_columns = ['id','code','pval_F','pval_S','pval_T','pval_F2','pval_S2','pval_T2','hits','similarity','score','align_length','overlap']

    with Pool() as pool, open(f'{file_name}_{domain}_{cluster}_comparison.csv','w', newline='') as csvfile:
        f_writer = csv.writer(csvfile)
        f_writer.writerow(data_columns)

        for ((key,code),p_value) in pool.imap_unordered(partial(calculatePvalue,remove_files=True),pdb_pairs,chunksize=50):
            if p_value != 'error':
                f_writer.writerow(['{0}_{1}_{4}_{2}_{3}_{5}'.format(*key),code]+[f'{n:.2e}' if isinstance(n,float) else str(n) for n in p_value])

def getDataset(heteromeric=True,filter_level=100,domain_type='CATH'):
    if not os.path.isfile(f'{"Heteromeric" if heteromeric else "Homomeric"}_complexes_{domain_type}_{filter_level}.csv'):
        print(f'Unable to find {"Heteromeric" if heteromeric else "Homomeric"} dataset, generating it now')
        makeDatasets(threshold=args.filter_level,use_identical_subunits=False,relabel=True,DEBUG_ignore_domains=False,domain_type=domain_type)
    print(f'Now loading dataset for {"Heteromeric" if heteromeric else "Homomeric"}_complexes_{domain_type}_{filter_level}.csv')
    return loadCSV(f'{"Heteromeric" if heteromeric else "Homomeric"}_complexes_{domain_type}_{filter_level}.csv')

def main(args):
    if args.filter_level not in VALID_CLUSTERS:
        print('Invalid filter level')
        return False

    heteromeric_data = getDataset(True,args.filter_level,args.domains)

    if args.stats:
        table, stats = shuffledInteractionDomains(heteromeric_data)
        print(f'Test statistic:\n\t{stats[0]:.3f}\np-value:\n\t{stats[1]:.3e}')
        print('Table:\n',table)

    if args.pullINT:
        pullXML(heteromeric_data['PDB_id'])

    if args.exec_mode == 'heteromeric':
        comparison_generator = heteromericInteractionRunner(heteromeric_data)
    elif args.exec_mode == 'precursor':
        homomeric_data = getDataset(False,args.filter_level,args.domains)
        if args.pullINT:
            pullXML(homomeric_data['PDB_id'])
        comparison_generator = precursorInteractionRunner(heteromeric_data,homomeric_data)
    else:
        print('Unknown comparison type')
        return False
        
    print('Running in parallel')
    paralleliseAlignment(comparison_generator,args.file_name,args.domains,args.filter_level)
    return

if __name__ == "__main__":
    os.environ['EMBOSS_ACDROOT'] = os.getcwd()
    parser = argparse.ArgumentParser(description = 'Domain comparison suite')

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-H", "--heteromeric", action="store_const",dest='exec_mode',const='heteromeric')
    group.add_argument("-P", "--precursor", action="store_const",dest='exec_mode',const='precursor')

    parser.add_argument('--filter', type=int,dest='filter_level',default=100)
    parser.add_argument('--file_name', type=str,dest='file_name',default='Data')
    parser.add_argument('--stats', action="store_true")
    parser.add_argument('--pullINT', action="store_true")
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--CATH", action="store_const",dest='domains',const='CATH')
    group.add_argument("--SCOP", action="store_const",dest='domains',const='SCOP')

    parser.set_defaults(exec_mode='heteromeric',domains='CATH')

    args = parser.parse_args()
    main(args)
