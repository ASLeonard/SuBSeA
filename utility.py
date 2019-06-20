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
