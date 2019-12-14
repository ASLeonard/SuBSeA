import pisa_XML, domains, SubSeA
import pytest
import os

def test_pull_XML(capsys):
    print('Starting .XML pulling and .INT generation test, stdout is enabled')
    with capsys.disabled():
        assert pisa_XML.pullXML(['4P69','1AI2','4WTO','1A1S','2PEY'])
    print('Cleaning out generated files')
    for ext_n in ['xml','int']:
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
    with capsys.disabled():
        SubSeA.calculatePvalue(('1A1S_A_A','1A1S_A_A','MUT'),WI_=False)
        print('With writing intermediate files')
        SubSeA.calculatePvalue(('1A1S_A_A','1A1S_A_A','MUT'),WI_=False)
    print('Test complete!\n')

def test_generate_datasets(capsys):
    #download xlxs
    #downdload clusters
    #strip and make
    pass



