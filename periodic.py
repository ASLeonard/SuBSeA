import pandas
from domains import readDomains, invertDomains

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

    result_rows = []

    for _, row in df.iterrows():
        pdb = row['PDB ID']

        domains = ALL_domains[pdb]

        ##empty domain, no point continuing
        if not domains:
            continue
        domain_edges = sorted(domains.keys())
        
        ordered_edges = []
        for r in row[key_interfaces].split(','):
            ordered_edges.extend(r.split('-'))
        
        ordered_edges = sorted(set(ordered_edges))
        
        interactions=[pair for pair, active in zip(row[key_interfaces].split(','),row[key_reduced].split(',')) if active == '1']


        for i,pair in enumerate(interactions):
            e1, e2 = pair.split('-')

            ##heteromeric interaction
            if e1 != e2:
                domains = ALL_domains[pdb]
                
                domain_edges = sorted(domains.keys())

                domain_overlap = 1 ##need to implement
                
                for edge in (e1, e2):
                    #comparables =
                    try:
                        edge = domain_edges[ordered_edges.index(edge)]
                    except:
                        print('bad length on', pdb)
                        print(ordered_edges, domain_edges)
                        continue
                    if tuple(domains[edge]) in inverted_homodimer_domains:
                        for comp in inverted_homodimer_domains[tuple(domains[edge])]:
                            yield (pdb + '_' + edge, '{}_{}'.format(*comp.split('_')))
                        #result_rows.append({'pdb':pdb, 'chain': edge, 'ref': comp, 'domain_overlap': 'full', 'heteromer_overlap': domain_overlap})
    return
        #for row in result_rows:
                    
                #find domain overlaps
                #run alignment on matches
                #store in nested heirarchy
                

    #return 1
def fnc(d):
    (het,hom) = pdb_combination
    args=(het[:4].upper(),het[5],hom[:4].upper(),hom[5])
    return '{}||{}'.format(args[0],args[2])
    
from multiprocessing import Pool
def qqq(gener):
    with Pool() as pool:
        results = pool.map(fnc,gener,50)
    return len(results)

def chainMap():
    chainmap = {}
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
        
if __name__ == "__main__":
    scrapePDBs(readHeteromers)
