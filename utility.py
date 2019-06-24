def formatFASTA(fname,out_name=None):
    if not out_name:
        last_index = fname.rfind('/') + 1
        out_name = fname[:last_index] + 'cleaned2_' + fname[last_index:]

    with open(fname,'r') as fasta_in, open(out_name,'w') as fasta_out:
        pdb, sequence = None, ''

        for line in fasta_in:
            if 'sequence' in line:
                pdb = '_'.join(line.split(':')[:2])
            elif 'secstr' in line:
                fasta_out.write(pdb+'\n'+sequence+'\n')
                pdb, sequence = None, ''
            elif pdb:
                sequence += line.rstrip()

from collections import defaultdict

def chainIt():
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
