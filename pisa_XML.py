#!/usr/bin/env python3                                                     
import sys
import os

from collections import defaultdict
import os.path

import xml.etree.ElementTree as ET

def splitXML(input_file='pisa.xml'):
    print('Spliting xml file')    
    file_string=''
    
    with open(input_file,'r') as file_:
        for line in file_:
            ##get the pdb id from the code tag
            if '<pdb_code>' in line:
                pdb_name = line.strip().split('>')[1].split('<')[0]
            ##magic string to indicate given XML depth of interest
            elif line[:2] == '  ' and line[:9] != '  <status':
                file_string += line[2:]
                ##end of entry, write to its own file
                if '</pdb_entry>' in line:
                    with open('XML/{}.xml'.format(pdb_name),'w') as out_file:
                        out_file.write(file_string)
                    file_string=''

def etree_to_dict(t):
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(('@' + k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
              d[t.tag]['#text'] = text
        else:
            d[t.tag] = text
    return d

def parseXML(xml_list = 'xml_list.txt'):
    ## consider all interfaces by default
    cssthresh = 0.0 

    with open(xml_list,'r') as file_:
        for pdb_entry in file_.readline().rstrip().split(', '):
            print('parsing entry '+pdb_entry)

            ##load tree in xml_format and convert recursively to dicts
            tree = ET.parse('XML/{}.xml'.format(pdb_entry))
            d = etree_to_dict(tree.getroot())

            if d['pdb_entry']['status'] == 'Ok':
                name = {}
                inter = {}
                if d['pdb_entry']['n_interfaces'] == '1':
                    d['pdb_entry']['interface'] = [d['pdb_entry']['interface']]
                for interface in d['pdb_entry']['interface']:
                    if float(interface['css']) < cssthresh or not all(molecule['class'] == 'Protein' for molecule in interface['molecule'][:2]):
                        continue

                    id_elements = ['{}_{}'.format(molecule['chain_id'],molecule['symop_no']) for molecule in interface['molecule'][:2]]
                    for idx,molecule in enumerate(interface['molecule'][:2]):
                        chain = molecule['chain_id']
                        if chain not in name:
                            name[chain]={}
                            inter[chain] = defaultdict(list)
                        symOP = molecule['symop_no']
                        residue = molecule['residues']['residue']
                        if isinstance(residue,dict):
                            residue=[residue]
                        for res in residue:
                            name[chain][int(res['seq_num'])] = res['name']
                            inter[chain][int(res['seq_num'])].append(id_elements[not idx]+' '+res['asa']+' '+res['bsa']+' '+id_elements[idx])                  

                ##write interactions to .int file
                with open('INT/{}.int'.format(pdb_entry),'w') as int_file:
                    for chain,residues in sorted(name.items()):
                        for res_seq,res_name in sorted(residues.items()):
                            int_file.write('{}\t{}\t{}\t'.format(chain,res_seq,res_name) + '\t'.join(inter[chain][res_seq])+'\n')
                        


if __name__ == "__main__":
    try:
        if len(sys.argv) == 1:
            print('Not enough arguments')
        elif sys.argv[1] == 'split':
            splitXML()
        elif '.xml' in sys.argv[1]:
            splitXML(sys.argv[1])
        elif sys.argv[1] == 'parse':
            parseXML()
        elif '.txt' in sys.argv[1]:
            parseXML(sys.argv[1])
        else:
            print('Unknown option')
    except Exception as e:
        print('Something went wrong with this',e)
