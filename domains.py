#!/usr/bin/env python3

import json
import requests
from collections import defaultdict, Counter


def pullDomains(pdb_id, style_key='CATH-B'):
     assert style_key == 'SCOP' or style_key == 'CATH-B', 'Unknown style_key: {}'.format(style_key)

     data=requests.get('http://www.ebi.ac.uk/pdbe/api/mappings/{}/{}'.format('scop' if style_key == 'SCOP' else 'cath_b', pdb_id),timeout=5).json()
     domain_arch=defaultdict(list)

     for class_id, superfam in data[pdb_id][style_key].items():
          if style_key == 'SCOP':
               class_id=superfam['superfamily']['sunid']
          for chain_element in superfam['mappings']:
               domain_arch[chain_element['chain_id']].append(class_id)
     return dict(domain_arch)

def writeDomains(pdb_list,class_type='SCOP',fname=None):
    if isinstance(pdb_list, str):
        with open(pdb_list,'r') as file_:
            pdb_list=file_.readline().rstrip().split(', ')
    ##pull domain information from EBI for all PDBs
    print('Pulling domains')
    domains={pdb.lower(): pullDomains(pdb.lower(),class_type) for pdb in pdb_list}

    ##write to json file
    print(f'writing to {fname}')
    with open(fname or f'domain_architectures_{class_type}.json', 'w') as file_out:
        file_out.write(json.dumps(domains))

def readDomains(file_name='SCOP'):
     with open(file_name if '.json' in file_name else f'domain_architectures_{file_name}.json') as f:
          return json.load(f)


def getUniqueHomodimerDomains(top_X=None):
     counted_domains=Counter([tuple(v) for vals in readDomains('HCath').values() for v in vals.values()])
     return counted_domains.most_common(top_X)


def domainMatch(data,c1,c2):
     if c1 not in data and c2 not in data:
          return 'NA'
     elif c1 not in data or c2 not in data:
          return 'none'
     elif data[c1] == data[c2]:
          return 'full'
     elif not set(data[c1]).isdisjoint(data[c2]):
          return 'partial'
     else:
          return 'none'

##invert domain data, into a dictionary with architecture key and pdb values
def invertDomains(domains,partials=False):
    inverted = defaultdict(list)

    for pdb_id, domain_info in domains.items():
        for chain, domains in domain_info.items():
            if isinstance(domains,list):
                domains = tuple(domains)
            if partials:
                for domain in set(domains):
                    inverted[domain].append(f'{pdb_id}_{chain}')
            else:
                inverted[domains].append(f'{pdb_id}_{chain}')

    return dict(inverted)
