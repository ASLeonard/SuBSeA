from SuBSeA import pisa_XML, domains, SubSeA, utility, pdb_visualisation
import os

def test_invert_domains(capsys):
    with capsys.disabled():
        print('Testing inverted domains, writing them first')
    domains.writeDomains(['4P69','1AI2','4WTO','1A1S','2PEY'],'CATH-B',fname='temporary.json')

    with capsys.disabled():
        print('Domains written, now reading in and inverting')
    domain_data = domains.readDomains('temporary.json')
    assert domains.invertDomains(domain_data)
    os.remove('temporary.json')

    with capsys.disabled():
        print('Test complete!\n')

#def test_needle():
#    print('Searching for needle exectuable')
    ## assert it exists in correct location

def test_SuBSeA(capsys):
    with capsys.disabled():
        print('Starting SuBSeA test, stdout is enabled')
        print('Getting files needed')

    pisa_XML.pullXML(['1A00','1BND'])
    for pdb, chains in (('1A00',('C','D')),('1BND',('A','B'))):
        for chain in chains:
            SubSeA.pullFASTA(pdb,chain)

    with capsys.disabled():
        print('Setting EMBOSS_ACDROOT variable to current directory')
        os.environ['EMBOSS_ACDROOT'] = os.getcwd()
        print('running calc')
        assert (SubSeA.calculatePvalue(('1A00_C_D','1A00_D_C','MUT'),WI_=False)[1] != 'error'), 'Didn\'t work'
        print('With writing intermediate files')
        assert (SubSeA.calculatePvalue(('1BND_A_B','1BND_B_A','MUT'),WI_=True)[1] != 'error'), 'Didn\'t work'
        print('calculations over, now cleaning files')
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

    with capsys.disabled():
        print('Test complete!\n')

def test_generate_datasets(capsys):
    with capsys.disabled():
        print('Making the dataset')
        utility.makeDatasets('',1,90,1,1,20)

    with capsys.disabled():
        print('Cleaning out files')
    os.remove('PeriodicTable.xlsx')
    os.remove('PDB_clusters_90.txt')

def test_heteromeric_grouping(capsys):
    with capsys.disabled():
        print('Including immunoglobulins')
    assert pdb_visualisation.groupHeteromericBSAs(90,False) is not None, 'Failed to group heteromeric interactions'

    with capsys.disabled():
        print('Excluding immunoglobulins')
    assert pdb_visualisation.groupHeteromericBSAs(90,True) is not None, 'Failed to group heteromeric interactions'
    os.remove('Heteromeric_complexes_90.csv')



