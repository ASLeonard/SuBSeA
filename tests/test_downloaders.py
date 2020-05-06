import sys, os
sys.path.append('SuBSeA')

from SuBSeA import pisa_XML, domains, binding_alignment, utility

def test_all_FASTA_download(capsys):
    with capsys.disabled():
        utility.downloadAllFASTA_GZ('')
    os.remove('minimal_all_fasta.txt')

def test_pull_XML(capsys):
    print('Starting .XML pulling and .INT generation test, stdout is enabled')
    with capsys.disabled():
        assert pisa_XML.pullXML(('4P69','1AI2','4WTO','1A1S','2PEY'))
    print('Cleaning out generated files')
    for temp_file in ('4P69','1AI2','4WTO','1A1S','2PEY'):
        os.remove(f'{temp_file}.int')
    print('Test complete!\n')

def test_pull_domains(capsys):
    print('Starting domain pulling test, stdout is enabled')
    with capsys.disabled():
        domains.writeDomains(['4P69','1AI2','4WTO','1A1S','2PEY'],'CATH-B',fname='temporary.json')
    print('Cleaning out generated files')
    os.remove('temporary.json')
    print('Test complete!\n')

def test_pull_fasta(capsys):
    print('Starting FASTA pulling test, stdout is enabled')
    with capsys.disabled():
        assert binding_alignment.pullFASTA('1A1S','A')
        assert binding_alignment.pullFASTA('4WTO','C')
    print('Cleaning out generated files')
    for temp_file in ('1A1S_A','4WTO_C'):
        os.remove(f'{temp_file}.fasta.txt')
    print('Test complete!\n')
