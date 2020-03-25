__all__ = ["DBlink"]
from dblink import DBlink

"""
PDB:
"""
PDBftp   = {'address'    : 'ftp.wwpdb.org',
            'structures' : '/pub/pdb/data/structures/divided/pdb/',
            'derived'    : '/pub/pdb/derived_data/index/',
            'resolution' : 'resolu.idx',
            'show'       : 'http://www.pdb.org/pdb/home/home.do'}
PDBrsync = {'port':   '33444',
            'address':'rsync.wwpdb.org::ftp/data/structures/divided/pdb/'}

"""
PDBeChem:
"""
PDBeChemftp = {'global':'ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif.tar.gz',
               'single':'ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif/',
               'show'  :'ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/'}

"""
GO:
"""
GOftp       = {'source':'ftp://ftp.geneontology.org/pub/go/ontology/go-simple.obo',
               'show'  : 'http://www.geneontology.org/'}
GOnamespace = {'biological_process':'B','cellular_component':'C','molecular_function':'M',
               'B':'biological_process','C':'cellular_component','M':'molecular_function'}

"""
TAXID:
"""
taxIDftp = {'global':'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip',
            'show'  :'http://www.ncbi.nlm.nih.gov/taxonomy'}

# """
# UNIPROT:
# """
Uniprotftp = {'swissprot':'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz',
               'trembl'   :'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz',
               'show'     :'http://www.uniprot.org/'}

"""
ENZYME:
"""
Enzymeftp = {'dat'  : 'ftp://ftp.expasy.org/databases/enzyme/enzyme.dat',
             'cls'  : 'ftp://ftp.expasy.org/databases/enzyme/enzclass.txt',
             'show' : 'http://enzyme.expasy.org/'}

"""
DRUGBANK:
"""
drugBankftp = {'show'    : 'http://www.drugbank.ca/',
               'main'    : 'http://www.drugbank.ca/system/downloads/current/drugbank.xml.zip',
               'targets' : 'http://www.drugbank.ca/system/downloads/current/all_target_ids_all.csv.zip'}

"""
SCOP:
"""
SCOPftp = {'show' : 'http://scop.mrc-lmb.cam.ac.uk/scop',
           'desc' : 'http://www.mrc-lmb.cam.ac.uk/agm/pre-scop/parseable/dir.des.scop.txt_MASTER.out',
           'rel'  : 'http://www.mrc-lmb.cam.ac.uk/agm/pre-scop/parseable/dir.cla.scop.txt_MASTER.out'}

"""
PDBTM:
"""
PDBTMftp = {'svn' : 'http://w2.enzim.hu/pdbtm/pdbtm',
            'show': 'http://pdbtm.enzim.hu/'}

"""
INCLUDES
"""
from PDBlink      import PDBlink
from PDBeChemlink import PDBeChemlink, PDBeChem
from GOlink       import GOlink, GOterm
from TaxIDlink    import TaxIDlink, TaxID
from Uniprotlink  import Uniprotlink, Uniprot
from Enzymelink   import Enzymelink, Enzyme
from DrugBanklink import DrugBanklink, Drug
from SCOPlink     import SCOPlink
from PDBTMlink    import PDBTMlink, TM

