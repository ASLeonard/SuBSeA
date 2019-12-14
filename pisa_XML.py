#!/usr/bin/env python3

import sys
import os

from collections import defaultdict

import xml.etree.ElementTree as ET # nosec

import urllib.request
import shutil
import argparse

BASE_PATH, XML_PATH, INT_PATH = '', '', ''

def checkINTExisting(df):
    missing_IDs = [row['PDB_id'] for _,row in clean_df.iterrows() if not os.path.exists(f'{BASE_PATH}{INT_PATH}/{row["PDB_id"].upper()}.int')]
    print(f'Mising {len(missing_IDs)} IDs')

def pullXML(pdb_code_file):
    print('Loading PDB files')
    print(BASE_PATH,XML_PATH,INT_PATH)
    
    pdbs = []

    if isinstance(pdb_code_file,str):
        with open(pdb_code_file) as file_in:
            for line in file_in:
                pdbs.extend(line.rstrip().split(', '))
    else:
        pdbs = pdb_code_file

    print(f'Loaded PDB files, there were {len(pdbs)} codes')

    URL_base = 'http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?'
    PER_CALL = 40

    print('Downloading XML from EBI')
    for slice_ in range(0,len(pdbs),PER_CALL):
        url=URL_base+','.join(map(str.lower,pdbs[slice_:slice_+PER_CALL]))
        with urllib.request.urlopen(url) as response, open(f'{BASE_PATH}{XML_PATH}XML_temp.xml', 'wb') as out_file: # nosec
            shutil.copyfileobj(response, out_file)
        splitXML(f'{BASE_PATH}{XML_PATH}XML_temp.xml')
        print(f'Split chunk {slice_} into XML_temp.xml')

    print('Cleaning temporary files')
    os.remove(f'{BASE_PATH}{XML_PATH}XML_temp.xml')

    print('Parsing individual files')
    parseXML(pdbs)
    print('Parsed all files')
    return True

def splitXML(input_file):
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
                    with open(f'{BASE_PATH}{XML_PATH}{pdb_name}.xml','w') as out_file:
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


##Converts XML file to INT file
def parseXML(xml_list):
    ## consider all interfaces by default
    cssthresh = 0.0

    if isinstance(xml_list, str):
        with open(xml_list,'r') as file_:
            xml_list=file_.readline().rstrip().split(', ')

    for pdb_entry in xml_list:
        pdb_entry=pdb_entry.upper()
        print(f'Parsing entry \'{pdb_entry}\'')

        ##load tree in xml_format and convert recursively to dicts
        try:
            tree = ET.parse(f'{BASE_PATH}{XML_PATH}{pdb_entry}.xml') # nosec
        except FileNotFoundError:
            print('Missing XML data on '+pdb_entry)
            continue
        d = etree_to_dict(tree.getroot())

        if d['pdb_entry']['status'] == 'Ok':
            name = {}
            inter = {}
            if d['pdb_entry']['n_interfaces'] == '1':
                d['pdb_entry']['interface'] = [d['pdb_entry']['interface']]
            for interface in d['pdb_entry']['interface']:
                if float(interface['css']) < cssthresh or not all(molecule['class'] == 'Protein' for molecule in interface['molecule'][:2]):
                    continue

                id_elements = [f"{molecule['chain_id']}_{molecule['symop_no']}" for molecule in interface['molecule'][:2]]
                for idx,molecule in enumerate(interface['molecule'][:2]):
                    chain = molecule['chain_id']
                    if chain not in name:
                        name[chain]={}
                        inter[chain] = defaultdict(list)
                    #symOP = molecule['symop_no']
                    residue = molecule['residues']['residue']
                    if isinstance(residue,dict):
                        residue=[residue]
                    for res in residue:
                        name[chain][int(res['seq_num'])] = res['name']
                        inter[chain][int(res['seq_num'])].append(id_elements[not idx]+' '+res['asa']+' '+res['bsa']+' '+id_elements[idx])

            ##write interactions to .int file
            with open(f'{BASE_PATH}{INT_PATH}{pdb_entry}.int','w') as int_file:
                print('writing to file')
                for chain,residues in sorted(name.items()):
                    for res_seq,res_name in sorted(residues.items()):
                        int_file.write(f'{chain}\t{res_seq}\t{res_name}\t' + '\t'.join(inter[chain][res_seq])+'\n')
        else:
            print(f'PDB entry {pdb_entry} status was not okay')

def main(args):

    global BASE_PATH
    BASE_PATH = args.file_path
    global XML_PATH
    XML_PATH = args.XML_path
    global INT_PATH
    INT_PATH = args.INT_path

    if args.file_format: ##default text file
        pullXML(args.file_name)
    else:
        pullXML(checkINTExisting(args.file_name))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'PDBePISA data scraper\nProduces .INT files for crystallographic interfaces.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--txt', action='store_true',dest='file_format')
    group.add_argument('--csv', action='store_false',dest='file_format')
    
    parser.add_argument('-N','--file_name', type=str,dest='file_name')
    parser.add_argument('-P','--file_path', type=str,dest='file_path')
    parser.add_argument('-X','--XML_path', type=str,dest='XML_path')
    parser.add_argument('-I','--INT_path', type=str,dest='INT_path')

    parser.set_defaults(file_format=True,file_name='data',file_path='../data/',XML_path='XML/',INT_path='INT/')

    args = parser.parse_args()
    
    main(args)
