#!/usr/bin/perl

use strict;
use warnings;

use autodie;
use Try::Tiny;

my $HOME_DIR = "/scratch/asl47/PDB";
my $WATER_EXEC = "/rscratch/asl47/water";
my $NEEDLE_EXEC = "/rscratch/asl47/needle";
my $FREESASA_EXEC = "/rscratch/asl47/freesasa";
my $TM_EXEC = "/rscratch/asl47/TMalign";

system("mkdir -p $HOME_DIR/results");

my ($pdb_raw, $HOMOMER, $KNOWN_BSA) = @ARGV;
my ($pdb_id, $BA_id) = split(/_/, $pdb_raw, 2);

my $subfolder = $HOMOMER eq 1 ? 'homodimer' : $HOMOMER;
my $results = "${HOME_DIR}/results/${subfolder}/${pdb_id}_${BA_id}.results";


print "Running DIMER analysis for $pdb_id with assembly $BA_id \n";
if (-f "$HOME_DIR/${pdb_id}.pdb${BA_id}") {
    #print "Existing pdb file\n"; 
}
else {
    #print "Downloading pdb file\n";
    system("wget -q \"https://files.rcsb.org/download/${pdb_id}.pdb${BA_id}.gz\" -P $HOME_DIR");
    #print "Downloaded fasta file\n";
    system("wget -q \"https://www.rcsb.org/pdb/download/downloadFastaFiles.do\?structureIdList=${pdb_id}&compressionType=uncompressed\" -O ${HOME_DIR}/${pdb_id}.fasta ");
    #print "Extracting pdb\n";	    
    system("gunzip -f $HOME_DIR/${pdb_id}.pdb${BA_id}.gz");

    #PROCESS PDB
    open (my $pdb_input, "<", "${HOME_DIR}/${pdb_id}.pdb${BA_id}");
    open (my $pdb_output, ">", "${HOME_DIR}/${pdb_id}_${BA_id}.pdb");

    {
        #local @ARGV = ($pdb_input);
        my $line_starts= "MODEL|ENDMDL|ATOM|TER|";

        while(<$pdb_input>){
            next unless /^MODEL.*|^ENDMDL.*|^ATOM.*|^TER.*/;
            print $pdb_output $_;
        }
        
       
    }
  
}

open (my $fasta, "<", "${HOME_DIR}/${pdb_id}.fasta");
try {
    my @chains;
    {
        local $/ = '>';
        while (<$fasta>) {
            s/^>//g; # strip out '>' from beginning
            s/>$//g; # and end of line
            next if !length($_);           # ignore empty lines

            my ($header_info) = /^(.*)\n/; # capture the header
            s/^(.*)\n//;                   # and strip it out
            my $chain_id=substr($header_info,5,1);
            #push @chains, $chain_id;
           
            s/\n//mg;

            open (my $fasta_chain, ">", "${HOME_DIR}/${pdb_id}.fasta${chain_id}");
            print $fasta_chain "$_";
            close $fasta_chain;
        }
    }
    open (my $pdb_innx, "<", "${HOME_DIR}/${pdb_id}_${BA_id}.pdb");
    while (<$pdb_innx>) {
        next unless /^TER.*/;       
        my @chain_id = split(' ',$_);
	my $chain_id_length=$chain_id[-1];
	if(scalar @chain_id ==3) 
	{
	$chain_id_length=substr($chain_id_length,1)+0;
	}
	else 
	{
	$chain_id_length =~ s/\D//g;
	#/-?\d+/
	}

	

        if(length($chain_id[2])==3 && $chain_id_length>20) {
            push @chains, substr($chain_id[3],0,1);
        }
        #print "CHN: " . @chain_id[3] . "\n";

    }
    
    system("echo ${HOMOMER} > ${results}"); 
    #print "LENGTH " . scalar @chains . "\n";
    if(scalar @chains != 2) {
        print "Bad length \n";
        goto OUTERLOOP; 
    }
    system("echo ${chains[0]} ${chains[1]} >> ${results}"); 
           
    if($HOMOMER eq "PDB_Homodimers") {
	system("${NEEDLE_EXEC} ${HOME_DIR}/${pdb_id}.fasta${chains[0]} ${HOME_DIR}/${pdb_id}.fasta${chains[1]} -gapopen 10 -gapextend .5 -nobrief -stdout -auto | grep Longest >> ${results}");
        system("${WATER_EXEC} ${HOME_DIR}/${pdb_id}.fasta${chains[0]} ${HOME_DIR}/${pdb_id}.fasta${chains[1]} -gapopen 10 -gapextend .5 -nobrief -stdout -auto | grep Longest >> ${results}");
        if($chains[0] ne $chains[1]) { 
            #print "Multi chain homomer\n";
            
            system("${FREESASA_EXEC} ${HOME_DIR}/${pdb_id}.pdb${BA_id} --join-models --chain-group=${chains[0]}${chains[1]}+${chains[0]}+${chains[1]} | grep Total | tail +2 >> ${results}");   

        }
        else {
            system("${FREESASA_EXEC} ${HOME_DIR}/${pdb_id}.pdb${BA_id} --join-models | grep Total >> ${results}");  
            system("${FREESASA_EXEC} ${HOME_DIR}/${pdb_id}.pdb${BA_id} --separate-models | grep Total >> ${results}"); 
        }
	system("~/Documents/PolyDev/duplication/dimer_split.py ${pdb_id} ${BA_id} ${chains[0]} ${chains[1]}");
 
        system("${TM_EXEC} ${HOME_DIR}/${pdb_id}_${chains[0]}.pdb ${HOME_DIR}/${pdb_id}_${chains[1]}.pdb -a | grep TM-score= | tail +3 >> ${results}");
    }
    else {
        #print "Heteromer\n";
	
	system("${NEEDLE_EXEC} ${HOME_DIR}/${pdb_id}.fasta${chains[0]} ${HOME_DIR}/${pdb_id}.fasta${chains[1]} -gapopen 10 -gapextend .5 -nobrief -stdout -auto | grep Longest >> ${results}");

        system("${WATER_EXEC} ${HOME_DIR}/${pdb_id}.fasta${chains[0]} ${HOME_DIR}/${pdb_id}.fasta${chains[1]} -gapopen 10 -gapextend .5 -nobrief -stdout -auto | grep Longest >> ${results}");
    
	if(!$KNOWN_BSA) {
        system("${FREESASA_EXEC} ${HOME_DIR}/${pdb_id}.pdb${BA_id} --chain-group=${chains[0]}${chains[1]}+${chains[0]}+${chains[1]} | grep Total | tail +2 >> ${results}");
	}   
        system("~/Documents/PolyDev/duplication/dimer_split.py ${pdb_id} ${BA_id} ${chains[0]} ${chains[1]}");
 
        system("${TM_EXEC} ${HOME_DIR}/${pdb_id}_${chains[0]}.pdb ${HOME_DIR}/${pdb_id}_${chains[1]}.pdb -a | grep TM-score= | tail +3 >> ${results}");
    }
} catch {
    warn "ID ${pdb_id} caught error: $_";
};

OUTERLOOP:
for my $file_d (grep {/fasta\S/} <${HOME_DIR}/*>) {
    #print "unlinking $file_d\n";
    unlink $file_d;
}

print "finished analysis of : $pdb_id\n";
