import pandas

def readHeteromers(file_path='~/Downloads/PeriodicTable.csv'):
    return pandas.read_csv(file_path)

def domainStrictness():
    pass

def scrapePDBs(df):
    with open('period_pdb_codes.txt', 'w') as file_out:
          file_out.write(', '.join(df['PDB ID']))

def trya(df):
    key='List of interface types (all identical subunits are given the same code)'

    N_SAMPLE_LIMIT = 500
    #row_iter = df.iterrows()

    result_rows=[]

    for _, row in df.iterrows():
        ##look at key
        pdb = row['PDB ID']
        interactions=[pair for pair, active in zip(row[5].split(','),row[8].split(',')) if active == '1']


        for i,pair in enumerate(interactions):
            e1, e2 = pair.split('-')

            ##heteromeric interaction
            if e1 != e2:
                domains = ALL_domains[pdb]
                domain_overlap = 1 ##need to implement
                
                for edge in (e1, e2):
                    comparables = inverted_homodimer_domains[domains[edge]]
                    for comp in comparables:
                        p_val = .05
                        result_rows.append({'pdb':pdb, 'chain': edge, 'ref': comp, 'domain_overlap': 'full', 'heteromer_overlap': domain_overlap})
                        
                
                #find domain overlaps
                #run alignment on matches
                #store in nested heirarchy
                1

    return 1

if __name__ == "__main__":
    scrapePDBs(readHeteromers)
