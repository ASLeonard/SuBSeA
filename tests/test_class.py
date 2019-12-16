import pisa_XML, domains, SubSeA
import pytest
import os

def test_pull_XML(capsys):
    print('Starting .XML pulling and .INT generation test, stdout is enabled')
    with capsys.disabled():
        assert pisa_XML.pullXML(('4P69','1AI2','4WTO','1A1S','2PEY'))
    print('Cleaning out generated files')
    for ext_n in ('xml','int'):
        for temp_file in ('4P69','1AI2','4WTO','1A1S','2PEY'):
            os.remove(f'{temp_file}.{ext_n}')
    print('Test complete!\n')

def test_pull_domains(capsys):
    print('Starting domain pulling test, stdout is enabled')
    with capsys.disabled():
        domains.writeDomains(['4P69','1AI2','4WTO','1A1S','2PEY'],'CATH-B',fname='temporary.json')
    print('Cleaning out generated files')
    os.remove('temporary.json')
    print('Test complete!\n')

def test_invert_domains():
    domains.writeDomains(['4P69','1AI2','4WTO','1A1S','2PEY'],'CATH-B',fname='temporary.json')
    domain_data = domains.readDomains('temporary.json')
    assert domains.invertDomains(domain_data)
    os.remove('temporary.json')
    print('Test complete!\n')

def test_pull_fasta(capsys):
    print('Starting FASTA pulling test, stdout is enabled')
    with capsys.disabled():
        assert SubSeA.pullFASTA('1A1S','A')
        assert SubSeA.pullFASTA('4WTO','C')
    print('Cleaning out generated files')
    for temp_file in ('1A1S_A','4WTO_C'):
        os.remove(f'{temp_file}.fasta.txt')
    print('Test complete!\n')

def test_needle():
    print('Searching for needle exectuable')
    ## assert it exists in correct location

def test_SuBSeA(capsys):
    print('Starting SuBSeA test, stdout is enabled')

    print('Getting files needed')
    pisa_XML.pullXML(['1A00','1BND'])
    for pdb, chains in (('1A00',('C','D')),('1BDN',('A','B'))):
        for chain in chains:
            SubSeA.pullFASTA(pdb,chain)
    
    with capsys.disabled():
        print('running calc')
        assert SubSeA.calculatePvalue(('1A00_C_D','1A00_D_C','MUT'),WI_=False)[1] != 'error', 'Didn\'t work'
        print('With writing intermediate files')
        SubSeA.calculatePvalue(('1BND_A_B','1BND_B_A','MUT'),WI_=True)[1] != 'error', 'Didn\'t work'

    ##clean up
    for pdb, chains in (('1A00',('C','D')),('1BND',('A','B'))):
        for chain in chains:
            os.remove(f'{pdb}_{chain}.fasta.txt')
    for ext_n in ('xml','int'):
        for temp_file in ('1A00','1BND'):
            os.remove(f'{temp_file}.{ext_n}')

    for ext_n in ('needle','pmatrix','pialign'):
        os.remove(f'1BND_A_1BND_B.{ext_n}')
    os.remove('1A00_C_1A00_D.needle')

    print('Test complete!\n')

def test_generate_datasets(capsys):
    #download xlxs
    #downdload clusters
    #strip and make
    pass



