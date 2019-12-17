from SubSeA import pullFASTA
from domains import pullDomains

from collections import defaultdict
import os
import json
import gzip
import shutil
import urllib.request
import ast


VALID_THRESHOLDS = {30,40,50,60,70,80,90,95,100}

##download file and then clean overall method
def downloadAllFASTA_GZ(fpath=''):
    print('Downloading compressed all seqres from PDB')
    urllib.request.urlretrieve('ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz', f'{fpath}all_fasta.txt.gz')
    
    print('Download successful, unzipping compressed archive now')
    with gzip.open(f'{fpath}all_fasta.txt.gz', 'rb') as f_in, open(f'{fpath}all_fasta.txt', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
        
    print('Cleaning up .txt.gz archive, now stripping out excess information')
    os.remove(f'{fpath}all_fasta.txt.gz')
    formatBulkFASTA(fpath,'all_fasta.txt')
    print(f'All FASTA sequences formatted into "{fpath}minimal_all_fasta.txt"')
    os.remove(f'{fpath}all_fasta.txt')
    

def formatBulkFASTA(fpath,fname,out_name=None):
    if not out_name:
        out_name = f'{fpath}minimal_{fname}'

    with open(fname,'r') as fasta_in, open(out_name,'w') as fasta_out:
        pdb, sequence = None, ''

        for line in fasta_in:
            if 'sequence' in line:
                pdb = '_'.join(line.split(':')[:2])
            elif 'secstr' in line:
                if '\n' in pdb or '\n' in sequence:
                    raise Exception('newline only line present in downloaded file, not sure how/why')
                if len(sequence) < 2:
                    print(pdb,sequence)
                    raise Exception(f'Something went wrong on {pdb}, {sequence}')
                fasta_out.write(f'{pdb}\n{sequence}\n')
                pdb, sequence = None, ''
            elif pdb:
                sequence += line.rstrip()

def loadCSV(fname):
    df = pandas.read_csv(fname,index_col=False)
    rr=[]
    for _,row in df.iterrows():

        interfaces = row['interfaces']
        if isinstance(interfaces,str):
            if '[' in interfaces:
                row['interfaces'] = ast.literal_eval(interfaces)
            else:
                row['interfaces'] = set(interfaces.split('-'))
        else:
            row['interfaces'] = ast.literal_eval(row['interfaces'].values[0])
        if '{' in row['domains']:
            row['domains'] = ast.literal_eval(row['domains'])
        else:
            row['domains'] = {dom.split(':')[0]:ast.literal_eval(dom.split(':')[1]) for dom in row['domains'].split(';') if len(dom)>1}
            

        row['BSAs'] = ast.literal_eval(row['BSAs'])
        #df.iloc[index] = row
        rr.append(row)
        
    return pandas.DataFrame(rr)

def writeCSV(df,fname):
    df.to_csv(fname,index=False,columns=['PDB_id','interfaces','domains','BSAs'])

def scrapePDBCodes(df):
    with open('periodic_pdb_codes.txt', 'w') as file_out:
        file_out.write(', '.join(df['PDB ID']))
    
from domains import readDomains, invertDomains

def invertCSVDomains(df, partials=False, homomeric=False):
    dom_dict = {}
    for _,row in df.iterrows():
        if row['domains'] is not None:
            if homomeric:
                dom_dict[row['PDB_id']] = {interaction.replace('-','_'):row['domains'][interaction[0]] for interaction in row['interfaces']}
            else:
                dom_dict[row['PDB_id']] = row['domains']
    return invertDomains(dom_dict,partials)
        
import pandas
from statistics import mean

def makeDatasets(fpath='',heteromeric=True, threshold=100):
    if not os.path.exists(f'{fpath}PeriodicTable.xlsx'):
        downloadPeriodicData(fpath)
    data = mergeSheets(fpath,heteromeric,1,1,0)
    post_filter = filterDataset(data,threshold)
    writeCSV(post_filter,('Heteromeric' if heteromeric else 'Homomeric')+f'_complexes_{threshold}.csv')

def downloadPeriodicData(fpath=''):
    print(f'Downloading periodic data from Science')
    urllib.request.urlretrieve('https://science.sciencemag.org/highwire/filestream/671215/' +
    'field_highwire_adjunct_files/3/aaa2245-Ahnert-SM-table-S2.xlsx', f'{fpath}PeriodicTable.xlsx')
    print(f'Download successful')

def mergeSheets(fpath='',heteromerics=True,use_identical_subunits=True,relabel=True,DEBUG_ignore_domains=True):

    DEX_interface = 'T' if heteromerics else 'I'

    interface_KEY = 'List of interface types (I- isologous homomeric; H - heterologous homomeric; T - heteromeric)'
    interface_LIST = 'List of interface types (all identical subunits are given the same code)' if use_identical_subunits else 'List of interfaces'
    interface_BSA= 'List of interface sizes (Angstroms^2)'

    assembly_SYM = ('Symmetry group (M - monomer; Dna - dihedral with heterologous interfaces;' +
    'Dns - dihedral with only isologous interfaces; Ts - tetrahedral with isologous interfaces;' +
    'Ta - tetrahedral with only heterologous interfaces; O* - 4 different octahedral topologies)')
    
    data_heteromers = pandas.read_excel(f'{fpath}PeriodicTable.xlsx',sheet_name=[0,2])

    try:
        domain_dict=readDomains('periodic2')
    except FileNotFoundError:
        print('given domain file doesn\'t exist, will start new one')
        domain_dict={}

    try:
        chain_map = chainMap() if relabel else {}
        chain_mapE = chainMap('Extra') if relabel else {}
        chain_map_full = {**chain_map,**chain_mapE}
    except FileNotFoundError:
        print('Chain relabelling doesn\'t exist, using blank')
        chain_map_full = {}
        
    new_rows = []

    for data in data_heteromers.values():
        for _, row in data.iterrows():
            PDB_code = row['PDB ID']

            if assembly_SYM in row and row[assembly_SYM] == 'M':
                continue
            
            
            if not pandas.isna(row[interface_KEY]) and any(type_ in row[interface_KEY] for type_ in ('T','I')):

                if not all(c.isupper() for pair in row[interface_KEY].split(',') for c in pair.split('-')):
                    print('werid chains',PDB_code,row[interface_KEY])
                    continue

                if PDB_code in chain_map_full:

                    new_interfaces = row[interface_LIST].lower()
                    for swaps in chain_map_full[PDB_code].items():
                        new_interfaces=new_interfaces.replace(swaps[0].lower(),swaps[1])
                    row[interface_LIST] = new_interfaces

                all_interfaces = zip(row[interface_LIST].split(','),row[interface_KEY].split(','))
                
                meaningful_interfaces = list({'-'.join(sorted(interface.split('-'))) for (interface,type_) in all_interfaces
                                            if (type_ == DEX_interface and (not heteromerics or interface[0] != interface[2]))})

                if not meaningful_interfaces:
                    continue

                BSAs = defaultdict(list)
                if not isinstance(row[interface_BSA],str):
                    row[interface_BSA] = str(row[interface_BSA])
                for K,V in zip(('-'.join(sorted(interface.split('-'))) for interface in row[interface_LIST].split(',')),row[interface_BSA].split(',')):
                    BSAs[K].append(float(V))

                BSA_av = {K:round(mean(V)) for K,V in BSAs.items()}

                if not DEBUG_ignore_domains and PDB_code not in domain_dict:
                    print('pulling for ',PDB_code)
                    domain_dict[PDB_code] = pullDomains(PDB_code)
                

                original_length = len(meaningful_interfaces)
                for index,interface in enumerate(meaningful_interfaces[:]):
                    if not all(chain in domain_dict[PDB_code] for chain in interface.split('-')):
                        del meaningful_interfaces[index + len(meaningful_interfaces) - original_length]
                if not meaningful_interfaces:
                    continue
                
                if DEBUG_ignore_domains:
                    domain_info = ''
                else:
                    domain_info_raw = ';'.join([chain+':{}'.format(tuple(domain_dict[PDB_code][chain])) if chain in domain_dict[PDB_code] 
                        else '' for chain in sorted({m for MI in meaningful_interfaces for m in MI.split('-')})])
                    domain_info = {dom.split(':')[0]:eval(dom.split(':')[1]) for dom in domain_info_raw.split(';') if len(dom)>1}
                
                new_rows.append({'PDB_id':row['PDB ID'], 'interfaces':meaningful_interfaces, 'domains':dsomain_info,
                    'BSAs': {K:BSA_av[K] for K in meaningful_interfaces}})

    if not DEBUG_ignore_domains:
        with open('domain_architectures_periodic2.json', 'w') as file_out:
            file_out.write(json.dumps(domain_dict))

    return pandas.DataFrame(new_rows)

## Download pre-made PDB clusters for protein subunits
def downloadClusters(threshold):
    assert threshold in VALID_THRESHOLDS, 'Invalid cluster threshold'
    print(f'Downloading PDB clustering @ {threshold}')
    urllib.request.urlretrieve(f'ftp://resources.rcsb.org/sequence/clusters/bc-{threshold}.out', f'PDB_clusters_{threshold}.txt')
    print(f'Download successful')

    
## Filter dataset at a specific redundancy level
def filterDataset(df,thresh,hom_mode=False):
    assert thresh in VALID_THRESHOLDS, f'Not implemented for threshold value: {thresh}'

    cluster_file_path = f'PDB_clusters_{thresh}.txt'

    if not os.path.exists(cluster_file_path):
        downloadClusters(thresh)
        
    with open(cluster_file_path) as cluster_file:
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
                print('Error, cannot be homomeric interaction but different clusters', pdb, interaction_pair)
            
            if cluster_indexes not in used_cluster_interactions:
                used_cluster_interactions.add(cluster_indexes)
                unique_interactions.append(interaction_pair)

        if unique_interactions:
            flat_chains = {C for pair in unique_interactions for C in pair.split('-')}
            new_df.append({'PDB_id':pdb,'interfaces':unique_interactions,'BSAs':{pair:row['BSAs'][pair] for pair in unique_interactions},'domains':{C:row['domains'][C] for C in flat_chains}})

    return pandas.DataFrame(new_df)

def chainMap(extra=None):
    chainmap = defaultdict(dict)
    with open('chain_map{}.txt'.format(extra or ''),'r') as file_:
        for line in file_:
            (pdb, file_chain, pdb_chain) = line.rstrip().split('\t')
            chainmap[pdb][file_chain] = pdb_chain
    return dict(chainmap)

def checkValidFASTA(file_name):
    return os.path.exists(file_name) and os.path.getsize(file_name)   
