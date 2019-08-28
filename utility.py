from SubSeA import pullFASTA
from domains import pullDomains
import json

def swapT(df,df2):
    rows = []
    for _,row in df.iterrows():
        row2 = df2.loc[df2['PDB_id']==row['PDB_id']]
        if row2 is not None:
            try:
                rows.append({'PDB_id':row['PDB_id'], 'interfaces':row['interfaces'],'domains':row['domains'],'BSAs':row2['BSAs'].values[0]})
            except:
                print(row2['BSAs'])
        else:
            print(row)
    return rows

def formatFASTA(fname,out_name=None):
    if not out_name:
        last_index = fname.rfind('/') + 1
        out_name = fname[:last_index] + 'xx_' + fname[last_index:]

    with open(fname,'r') as fasta_in, open(out_name,'w') as fasta_out:
        pdb, sequence = None, ''

        for line in fasta_in:
            if 'sequence' in line:
                pdb = '_'.join(line.split(':')[:2])
            elif 'secstr' in line:
                if '\n' in pdb or '\n' in sequence:
                    raise Excepion('f')
                if len(sequence) < 2:
                    print(pdb,sequence)
                    raise Exception('Something went wrong!')
                fasta_out.write(pdb+'\n'+sequence+'\n')
                pdb, sequence = None, ''
            elif pdb:
                sequence += line.rstrip()


def loadCSV(fname):
    df = pandas.read_csv(fname,index_col=False)
    rr=[]
    for index,row in df.iterrows():
        #row['BSAs']

        interfaces = row['interfaces']
        if isinstance(interfaces,str):
            if '[' in interfaces:
                row['interfaces'] = eval(interfaces)
            else:
                row['interfaces'] = set(interfaces.split('-'))
        else:
            row['interfaces'] = eval(row['interfaces'].values[0])
        if True or ';' in row['domains']:
            row['domains'] = {dom.split(':')[0]:eval(dom.split(':')[1]) for dom in row['domains'].split(';') if len(dom)>1}
        else:
            row['domains'] = eval(row['domains'])

        row['BSAs'] = eval(row['BSAs'])
        #df.iloc[index] = row
        rr.append(row)
        
    return pandas.DataFrame(rr)

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

def getBadIDs(dx=None,raw_table=None):
    if raw_table is None:
        raw_table  = pandas.read_csv('Periodic_heteromers3.csv',index_col=False)
    if not dx:
        ps = list(raw_table['PDB_id'])
        lx= [] 
        for line in open('/scratch/asl47/PDB/FASTA/cleaned2_all_fasta.txt','r'):
            if '>' in line and line.rstrip()[1:5].lower() in ps:
                lx.append(line.rstrip()[1:])

        dx = defaultdict(list)
        for ix in lx:
            dx[ix[:4].lower()].append(ix[-1])
        return dx
    
    bad_s=[]

    #raw_table  = pandas.read_csv('Periodic_heteromers3.csv',index_col=False)
    for _,row in raw_table.iterrows():
        #row['interfaces'] = eval(row['interfaces'].values[0])
        #row['interfaces'] = {i for i in row['interfaces'][2:-2].split('\', \'')}
        interfaces = {chain for interaction in row['interfaces'] for chain in interaction.split('-')}
        for inx in interfaces:
            if inx not in dx[row['PDB_id']]:
                bad_s.append('{}_{}'.format(row['PDB_id'],inx))
    return bad_s

def nextWave(bads):
    super_bad = []
    for bad in bads:
        if not fullFASTA('/scratch/asl47/PDB/FASTA/{}.fasta.txt'.format(bad.upper())):
            pullFASTA(*bad.split('_'))
            if not fullFASTA('/scratch/asl47/PDB/FASTA/{}.fasta.txt'.format(bad.upper())):
                super_bad.append(bad)
    return super_bad
                
        

def mergeSheets(heteromerics=True):

    DEX_interface = 'T' if heteromerics else 'I'

    interface_KEY = 'List of interface types (I- isologous homomeric; H - heterologous homomeric; T - heteromeric)'
    interface_LIST = 'List of interface types (all identical subunits are given the same code)'
    interface_BSA= 'List of interface sizes (Angstroms^2)'

    assembly_SYM = 'Symmetry group (M - monomer; Dna - dihedral with heterologous interfaces; Dns - dihedral with only isologous interfaces; Ts - tetrahedral with isologous interfaces; Ta - tetrahedral with only heterologous interfaces; O* - 4 different octahedral topologies)'

    redundant = 'Included in non-redundant set'
    
    data_heteromers = pandas.read_excel('~/Downloads/PeriodicTable.xlsx',sheet_name=[2])

    domain_dict=readDomains('periodic3')
    chain_map = chainMap()
    chain_mapE = chainMap('Extra')
    new_rows = []

    for data in data_heteromers.values():
        for i, row in data.iterrows():
            if redundant in row and row[redundant] == 0:
                continue
            PDB_code = row['PDB ID']
            
            if assembly_SYM in row and row[assembly_SYM] == 'M':
                continue
            
            
            if not pandas.isna(row[interface_KEY]) and any(type_ in row[interface_KEY] for type_ in ('T','I')):

                if PDB_code in chain_map:
                    for swaps in chain_map[PDB_code].items():
                        row[interface_LIST] = row[interface_LIST].replace(*swaps)
                        
                if PDB_code in chain_mapE:
                    for swaps in chain_mapE[PDB_code].items():
                        row[interface_LIST] = row[interface_LIST].replace(*swaps)
                        
                all_interfaces = zip(row[interface_LIST].split(','),row[interface_KEY].split(','))   
                
                meaningful_interfaces = list({'-'.join(sorted(interface.split('-'))) for (interface,type_) in all_interfaces if (type_ == DEX_interface and (not heteromerics or interface[0] != interface[2]))})

                if not meaningful_interfaces:
                    continue

                BSAs = defaultdict(list)
                if not isinstance(row[interface_BSA],str):
                    row[interface_BSA] = str(row[interface_BSA])
                for K,V in zip(('-'.join(sorted(interface.split('-'))) for interface in row[interface_LIST].split(',')),row[interface_BSA].split(',')):
                    BSAs[K].append(float(V))

                BSA_av = {K:round(mean(V)) for K,V in BSAs.items()}

                if PDB_code not in domain_dict:
                    print('pulling for ',PDB_code)
                    domain_dict[PDB_code] = pullDomains(PDB_code)

                original_length = len(meaningful_interfaces)
                for index,interface in enumerate(meaningful_interfaces[:]):
                    if not all(chain in domain_dict[PDB_code] for chain in interface.split('-')):
                        del meaningful_interfaces[index + len(meaningful_interfaces) - original_length]
                if not meaningful_interfaces:
                    continue
                
                domain_info = ';'.join([chain+':{}'.format(tuple(domain_dict[PDB_code][chain])) if chain in domain_dict[PDB_code] else '' for chain in sorted({m for MI in meaningful_interfaces for m in MI.split('-')})])
                
                new_rows.append({'PDB_id':row['PDB ID'], 'interfaces':meaningful_interfaces, 'domains':domain_info, 'BSAs': str({K:BSA_av[K] for K in meaningful_interfaces})})

                
    with open('domain_architectures_periodic3.json', 'w') as file_out:
        file_out.write(json.dumps(domain_dict))

    return pandas.DataFrame(new_rows)


def chainMap(extra=None):
    chainmap = defaultdict(dict)
    with open('chain_map{}.txt'.format(extra or ''),'r') as file_:
        for line in file_:
            (pdb, file_chain, pdb_chain) = line.rstrip().split('\t')
            chainmap[pdb][file_chain] = pdb_chain
    return dict(chainmap)

from collections import defaultdict

def fullFASTA(file_name):
    return os.path.exists(file_name) and os.path.getsize(file_name)


def chainFull():
    chainmap = {}
    fc = open('chain_map.txt')
    for line in fc:
        l = line.strip().split('\t')
        if l[0].upper() not in chainmap:
            chainmap[l[0].upper()] = {}
        chainmap[l[0].upper()][l[1]] = l[2]

    chainmapset = set(chainmap.keys())

    chainset = defaultdict(list)
    fc2 = open('pdb_chains.txt')
    for line in fc2:
        l = line.strip().split('\t')
        chainset[l[0].upper()].append(set(l[1].split()))
    #return chainset
        
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
                pfam[i][j] = ';'.join(sorted(list(pfaml[i][j])))
    return chainmap,chainmapset,pfaml,pfam


def fixEm(pchains,fixables,df):
    rev_map=defaultdict(dict)
    for fix in fixables:
        inters = eval(df.loc[df['PDB_id']==fix[:4]]['interfaces'].values[0])
        #print(inters)
        #return
        pch = pchains[fix[:4].upper()]
        pch = [i for j in  pch for i in j]
        inters2= [i for j in inters for i in j.split('-')]
        #print(pch,inters2)
        if len(inters2) == 2 and len(pch)==2:
            #return fix
            #dif1 = set(inters2)-set(pch)
            #dif2 = set(pch) - set(inters2)
            #if len(dif1) == 1 and len(dif2) ==1:
            #    return fix
            for i in range(2):
                rev_map[fix[:4]][inters2[i]] = pch[i]
            #rev_map[fix[:4]][inters2]=list(pch[0])[0]
            #rev_map[fix[:4]][inters[inters.index('-')+1]]=list(pch[1])[0]
    return rev_map

def wrFE(rm):
    #rev_map[fix[:4]][inters[inters.index('-')-1]]=list(pch[0])[0]
    tw= ''
    for pdb, maps in rm.items():
        
        for A,B in maps.items():
            tw+='\t'.join((pdb,A,B))+'\n'
    return tw
