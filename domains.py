#!/usr/bin/env python3
import sys
import json
import requests
from collections import defaultdict, Counter


def pullDomains(pdb_id, style_key='CATH-B'):
     assert style_key == 'SCOP' or style_key == 'CATH-B', 'Unknown style_key: {}'.format(style_key)

     data=requests.get('http://www.ebi.ac.uk/pdbe/api/mappings/{}/{}'.format('scop' if style_key == 'SCOP' else 'cath_b', pdb_id)).json()
     domain_arch=defaultdict(list)

     for class_id, superfam in data[pdb_id][style_key].items():
          if style_key == 'SCOP':
               class_id=superfam['superfamily']['sunid']
          for chain_element in superfam['mappings']:
               domain_arch[chain_element['chain_id']].append(class_id)
     return dict(domain_arch)

def writeDomains(pdb_list,class_type='SCOP'):
    if isinstance(pdb_list, str):
        with open(pdb_list,'r') as file_:
            pdb_list=file_.readline().rstrip().split(', ')
    ##pull domain information from EBI for all PDBs
    print('Pulling domains')
    domains={pdb: pullDomains(pdb,class_type) for pdb in pdb_list}

    ##write to json file
    with open('domain_architectures_{}.json'.format(class_type), 'w') as file_out:
        file_out.write(json.dumps(domains))

def readDomains(file_name='SCOP'):
     with open('domain_architectures_{}.json'.format(file_name)) as f:
          return json.load(f)


def getUniqueHomodimerDomains(top_X=None):
     if top_X is None:
          return {tuple(v) for vals in readDomains('HCath').values() for v in vals.values()}
     else:
          counted_domains=Counter([tuple(v) for vals in readDomains('HCath').values() for v in vals.values()])
          return {mc[0] for mc in counted_domains.most_common(top_X)}


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
def invertDomains(domains,ids=None):
     inverted=defaultdict(list)

     for k,v in domains.items():
          if not v:
               continue
          #elif k not in ids:
          #     continue
          #elif ids and len(set(tuple(v[chain]) for chain in ids[k])):
          #     continue
          arch=tuple(list(v.values())[0])
          inverted[arch].append('{}_{}'.format(k,list(v.keys())[0]))
     return dict(inverted)
if __name__ == "__main__":
    try:
        if len(sys.argv) == 1:
            print('Not enough arguments')
        elif '.txt' in sys.argv[1]:
            print('About to pull domains')
            writeDomains(sys.argv[1],'CATH-B')
    except Exception as e:
        print('err',e)
