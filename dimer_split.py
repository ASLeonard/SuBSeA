#!/usr/bin/python3

import sys
import os


def analyseMany(f_name,dir_out):
    if 'PDB' in f_name:
        for line in [l.rstrip().split(', ') for l in open('Inputs/{}.txt'.format(f_name))][0]:
            os.system('/u/fs1/asl47/Documents/PolyDev/duplication/dimer_analysis.pl {}_1 {} 0'.format(line,dir_out))
    else:
        for line in open('Inputs/{}.txt'.format(f_name)):
            (pdb_id,BA_id) =str(line.split('\t')[0]).split('_')
            os.system('/u/fs1/asl47/Documents/PolyDev/duplication/dimer_analysis.pl {}_{} {} 0'.format(pdb_id,BA_id,dir_out))
        
        

def splitChainFiles(pdb_in,pdb_ID, ch1, ch2):
    with open('{}_{}.pdb'.format(pdb_in,ch1), 'w') as ch1_out:
        with open('{}_{}.pdb'.format(pdb_in,ch2), 'w') as ch2_out:
            for line in open('{}.pdb{}'.format(pdb_in,pdb_ID)):
                if 'ATOM' in line:
                    if ch1 in line.split()[:-1]:
                        ch1_out.write(line)
                    elif ch1!=ch2 and ch2 in line.split()[:-1]:
                        ch2_out.write(line)

def main():
    if len(sys.argv)==5:
        splitChainFiles(*sys.argv[1:5])
    elif len(sys.argv)==2:
        analyseMany(sys.argv[1],sys.argv[1])
        print('Finished')
    else:
        print("No args")
            
if __name__ == "__main__":
    os.chdir('/scratch/asl47/PDB')
    main()
