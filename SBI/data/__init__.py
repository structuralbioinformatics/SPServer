'''
#
# PROTEIN
#
'''

'''
CODING TRANSFORMATION
http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/
'''

aminoacids3to1 = {
    'ALA': 'A', 'AZT': 'A', 'CHA': 'A', 'HPH': 'A', 'NAL': 'A', 'AIB': 'A', 'BAL': 'A',
    'DHA': 'A', 'BB9': 'A', 'ALM': 'A', 'AYA': 'A', 'BNN': 'A', 'CHG': 'A', 'CSD': 'A',
    'DAL': 'A', 'DNP': 'A', 'FLA': 'A', 'HAC': 'A', 'MAA': 'A', 'PRR': 'A', 'TIH': 'A',
    'TPQ': 'A', 'BB9': 'A',
    'ARG': 'R', 'ORN': 'R', 'ACL': 'R', 'ARM': 'R', 'AGM': 'R', 'HAR': 'R', 'HMR': 'R',
    'DAR': 'R',
    'ASN': 'N', 'MEN': 'N',
    'ASP': 'D', 'ASZ': 'D', '2AS': 'D', 'ASA': 'D', 'ASB': 'D', 'ASK': 'D', 'ASL': 'D',
    'ASQ': 'D', 'BHD': 'D', 'DAS': 'D', 'DSP': 'D',
    'ASX': 'B',
    'CYS': 'C', 'CYD': 'C', 'CYO': 'C', 'HCY': 'C', 'CSX': 'C', 'SMC': 'C', 'BCS': 'C',
    'BUC': 'C', 'C5C': 'C', 'C6C': 'C', 'CCS': 'C', 'CEA': 'C', 'CME': 'C', 'CSO': 'C',
    'CSP': 'C', 'CSS': 'C', 'CSW': 'C', 'CY1': 'C', 'CY3': 'C', 'CYG': 'C', 'CYM': 'C',
    'CYQ': 'C', 'DCY': 'C', 'OCS': 'C', 'SOC': 'C', 'EFC': 'C', 'PR3': 'C', 'SCH': 'C',
    'SCS': 'C', 'SCY': 'C', 'SHC': 'C', 'PEC': 'C',
    'GLN': 'Q', 'DGN': 'Q',
    'GLU': 'E', 'GLA': 'E', 'GLZ': 'E', 'PCA': 'E', '5HP': 'E', 'CGU': 'E', 'DGL': 'E',
    'GGL': 'E', 'GMA': 'E',
    'GLX': 'Z',
    'GLY': 'G', 'GL3': 'G', 'GLZ': 'G', 'GSC': 'G', 'SAR': 'G', 'MPQ': 'G', 'NMC': 'G',
    'MSA': 'G', 'DBU': 'G',
    'HIS': 'H', 'HSE': 'H', 'HSD': 'H', 'HI0': 'H', 'HIP': 'H', 'HID': 'H', 'HIE': 'H',
    '3AH': 'H', 'MHS': 'H', 'DHI': 'H', 'HIC': 'H', 'NEP': 'H', 'NEM': 'H',
    'ILE': 'I', 'IIL': 'I', 'DIL': 'I',
    'LEU': 'L', 'NLE': 'L', 'LOV': 'L', 'NLN': 'L', 'NLP': 'L', 'MLE': 'L', 'BUG': 'L',
    'CLE': 'L', 'DLE': 'L', 'MLU': 'L',
    'LYS': 'K', 'LYZ': 'K', 'ALY': 'K', 'TRG': 'K', 'SHR': 'K', 'LYM': 'K', 'LLY': 'K',
    'KCX': 'K', 'LLP': 'K', 'DLY': 'K', 'DM0': 'K',
    'MET': 'M', 'MSE': 'M', 'CXM': 'M', 'FME': 'M', 'OMT': 'M',
    'PHE': 'F', 'DAH': 'F', 'HPQ': 'F', 'DPN': 'F', 'PHI': 'F', 'PHL': 'F',
    'PRO': 'P', 'HYP': 'P', 'DPR': 'P', 'ECQ': 'P', 'POM': 'P', 'H5M': 'P',
    'SER': 'S', 'HSE': 'S', 'STA': 'S', 'SVA': 'S', 'SAC': 'S', 'SEL': 'S', 'SEP': 'S',
    'SET': 'S', 'OAS': 'S', 'DSN': 'S', 'MIS': 'S',
    'THR': 'T', 'PTH': 'T', 'ALO': 'T', 'TPO': 'T', 'BMT': 'T', 'DTH': 'T', 'CTH': 'T',
    'TRP': 'W', 'TPL': 'W', 'TRO': 'W', 'DTR': 'W', 'HTR': 'W', 'LTR': 'W',
    'TYR': 'Y', 'TYQ': 'Y', 'TYS': 'Y', 'TYY': 'Y', 'TYB': 'Y', 'STY': 'Y', 'PTR': 'Y',
    'PAQ': 'Y', 'DTY': 'Y', 'IYR': 'Y', 'GHP': 'Y', 'D3P': 'Y', 'D4P': 'Y', 'OMZ': 'Y',
    'OMY': 'Y',
    'VAL': 'V', 'NVA': 'V', 'DVA': 'V', 'DIV': 'V', 'MVA': 'V',
    'SEC': 'U',
    'PYL': 'O',
    'XLE': 'J',
    'ACE': 'X', '3FG': 'X', 'UNK': 'X'
}

aminoacids1to3 = dict([[v, k] for k, v in aminoacids3to1.items()])
aminoacids1to3['A'] = 'ALA'
aminoacids1to3['N'] = 'ASN'
aminoacids1to3['R'] = 'ARG'
aminoacids1to3['D'] = 'ASP'
aminoacids1to3['C'] = 'CYS'
aminoacids1to3['Q'] = 'GLN'
aminoacids1to3['E'] = 'GLU'
aminoacids1to3['G'] = 'GLY'
aminoacids1to3['H'] = 'HIS'
aminoacids1to3['I'] = 'ILE'
aminoacids1to3['J'] = 'XLE'
aminoacids1to3['L'] = 'LEU'
aminoacids1to3['K'] = 'LYS'
aminoacids1to3['M'] = 'MET'
aminoacids1to3['F'] = 'PHE'
aminoacids1to3['O'] = 'PYL'
aminoacids1to3['P'] = 'PRO'
aminoacids1to3['S'] = 'SER'
aminoacids1to3['T'] = 'THR'
aminoacids1to3['U'] = 'SEC'
aminoacids1to3['W'] = 'TRP'
aminoacids1to3['Y'] = 'TYR'
aminoacids1to3['V'] = 'VAL'

'''
REGULAR AMINOACIDS IDENTIFICATION
'''
aminoacids_main3 = set(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU',
                        'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'SEC', 'PYL', 'XLE'])

aminoacids_main1 = set(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])

'''
PROPERTIES
'''
aminoacids_surface = {
    'A': 115, 'C': 149, 'D': 170, 'E': 207, 'F': 230, 'G': 86,  'H': 206,
    'I': 187, 'K': 222, 'L': 192, 'M': 210, 'N': 184, 'P': 140, 'Q': 208,
    'R': 263, 'S': 140, 'T': 164, 'V': 161, 'W': 269, 'Y': 257,
}

aminoacids_polarity_boolean = {
    'A': False, 'C': False, 'D': True,  'E': True,  'F': False, 'G': False, 'H': True,
    'I': False, 'K': True,  'L': False, 'M': False, 'N': True,  'P': False, 'Q': True,
    'R': True,  'S': True,  'T': True,  'V': False, 'W': False, 'Y': True
}

'''
#
# DNA/RNA
#
'''
nucleic_main = set(['DA', 'A', 'DC', 'C', 'DG', 'G', 'DI', 'I', 'DT', 'T', 'DU', 'U'])
nucleic2to1  = {
    'DA': 'A', 'A': 'A', 'A44': 'A', '6HA': 'A', 'APN': 'A',
    'DC': 'C', 'C': 'C', '5CM': 'C', 'MCY': 'C', 'OMC': 'C', 'C43': 'C', '6HC': 'C',
    'CPN': 'C',
    'DG': 'G', 'G': 'G', 'OMG': 'G', 'G48': 'G', '6HG': 'G', 'GPN': 'G',
    'DI': 'I', 'I': 'I',
    'DT': 'T', 'T': 'T', '6HT': 'T', 'TPN': 'T',
    'DU': 'U', 'U': 'U', 'U36': 'U', '5IU': 'U',
    'N': 'N'
}

'''
CRYSTALOGRAPHIC METHODS
'''
crystal_method_has_resolution = set(['X-RAY DIFFRACTION', 'ELECTRON MICROSCOPY', 'NEUTRON DIFFRACTION',
                                     'FIBER DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY'])

crystal_method_not_resolution = set(['SOLUTION NMR', 'POWDER DIFFRACTION', 'SOLUTION SCATTERING', 'SOLID-STATE NMR',
                                     'INFRARED SPECTROSCOPY', 'FLUORESCENCE TRANSFER'])

crystal_method                = crystal_method_has_resolution.union(crystal_method_not_resolution)


'''
Elements
'''


class Element(object):
    def __init__(self, number, symbol, name):
        self.number = number
        self.symbol = symbol
        self.name   = name

element_dic = {
    'H':  Element(  1, 'H',  'Hydrogen'),    'He': Element(  2, 'He', 'Helium'),
    'Li': Element(  3, 'Li', 'Lithium'),     'Be': Element(  4, 'Be', 'Beryllium'),
    'B':  Element(  5, 'B',  'Boron'),       'C':  Element(  6, 'C',  'Carbon'),
    'N':  Element(  7, 'N',  'Nitrogen'),    'O':  Element(  8, 'O',  'Oxygen'),
    'F':  Element(  9, 'F',  'Fluorine'),    'Ne': Element( 10, 'Ne', 'Neon'),
    'Na': Element( 11, 'Na', 'Sodium'),      'Mg': Element( 12, 'Mg', 'Magnesium'),
    'Al': Element( 13, 'Al', 'Aluminium'),   'Si': Element( 14, 'Si', 'Silicon'),
    'P':  Element( 15, 'P',  'Phosphorus'),  'S':  Element( 16, 'S',  'Sulfur'),
    'Cl': Element( 17, 'Cl', 'Chlorine'),    'Ar': Element( 18, 'Ar', 'Argon'),
    'K':  Element( 19, 'K',  'Potassium'),   'Ca': Element( 20, 'Ca', 'Calcium'),
    'Sc': Element( 21, 'Sc', 'Scandium'),    'Ti': Element( 22, 'Ti', 'Titanium'),
    'V':  Element( 23, 'V',  'Vanadium'),    'Cr': Element( 24, 'Cr', 'Chromium'),
    'Mn': Element( 25, 'Mn', 'Manganese'),   'Fe': Element( 26, 'Fe', 'Iron'),
    'Co': Element( 27, 'Co', 'Cobalt'),      'Ni': Element( 28, 'Ni', 'Nickel'),
    'Cu': Element( 29, 'Cu', 'Copper'),      'Zn': Element( 30, 'Zn', 'Zinc'),
    'Ga': Element( 31, 'Ga', 'Gallium'),     'Ge': Element( 32, 'Ge', 'Germanium'),
    'As': Element( 33, 'As', 'Arsenic'),     'Se': Element( 34, 'Se', 'Selenium'),
    'Br': Element( 35, 'Br', 'Bromine'),     'Kr': Element( 36, 'Kr', 'Krypton'),
    'Rb': Element( 37, 'Rb', 'Rubidium'),    'Sr': Element( 38, 'Sr', 'Strontium'),
    'Y':  Element( 39, 'Y',  'Yttrium'),     'Zr': Element( 40, 'Zr', 'Zirconium'),
    'Nb': Element( 41, 'Nb', 'Niobium'),     'Mo': Element( 42, 'Mo', 'Molybdenum'),
    'Tc': Element( 43, 'Tc', 'Technetium'),  'Ru': Element( 44, 'Ru', 'Ruthenium'),
    'Rh': Element( 45, 'Rh', 'Rhodium'),     'Pd': Element( 46, 'Pd', 'Palladium'),
    'Ag': Element( 47, 'Ag', 'Silver'),      'Cd': Element( 48, 'Cd', 'Cadmium'),
    'In': Element( 49, 'In', 'Indium'),      'Sn': Element( 50, 'Sn', 'Tin'),
    'Sb': Element( 51, 'Sb', 'Antimony'),    'Te': Element( 52, 'Te', 'Tellurium'),
    'I':  Element( 53, 'I',  'Iodine'),      'Xe': Element( 54, 'Xe', 'Xenon'),
    'Cs': Element( 55, 'Cs', 'Caesium'),     'Ba': Element( 56, 'Ba', 'Barium'),
    'La': Element( 57, 'La', 'Lanthanum'),   'Ce': Element( 58, 'Ce', 'Cerium'),
    'Pr': Element( 59, 'Pr', 'Praseodymium'), 'Nd': Element( 60, 'Nd', 'Neodymium'),
    'Pm': Element( 61, 'Pm', 'Promethium'),  'Sm': Element( 62, 'Sm', 'Samarium'),
    'Eu': Element( 63, 'Eu', 'Europium'),    'Gd': Element( 64, 'Gd', 'Gadolinium'),
    'Tb': Element( 65, 'Tb', 'Terbium'),     'Dy': Element( 66, 'Dy', 'Dysprosium'),
    'Ho': Element( 67, 'Ho', 'Holmium'),     'Er': Element( 68, 'Er', 'Erbium'),
    'Tm': Element( 69, 'Tm', 'Thulium'),     'Yb': Element( 70, 'Yb', 'Ytterbium'),
    'Lu': Element( 71, 'Lu', 'Lutetium'),    'Hf': Element( 72, 'Hf', 'Hafnium'),
    'Ta': Element( 73, 'Ta', 'Tantalum'),    'W':  Element( 74, 'W', 'Tungsten'),
    'Re': Element( 75, 'Re', 'Rhenium'),     'Os': Element( 76, 'Os', 'Osmium'),
    'Ir': Element( 77, 'Ir', 'Iridium'),     'Pt': Element( 78, 'Pt', 'Platinum'),
    'Au': Element( 79, 'Au', 'Gold'),        'Hg': Element( 80, 'Hg', 'Mercury'),
    'Tl': Element( 81, 'Tl', 'Thallium'),    'Pb': Element( 82, 'Pb', 'Lead'),
    'Bi': Element( 83, 'Bi', 'Bismuth'),     'Po': Element( 84, 'Po', 'Polonium'),
    'At': Element( 85, 'At', 'Astatine'),    'Rn': Element( 86, 'Rn', 'Radon'),
    'Fr': Element( 87, 'Fr', 'Francium'),    'Ra': Element( 88, 'Ra', 'Radium'),
    'Ac': Element( 89, 'Ac', 'Actinium'),    'Th': Element( 90, 'Th', 'Thorium'),
    'Pa': Element( 91, 'Pa', 'Protactinium'), 'U':  Element( 92, 'U', 'Uranium'),
    'Np': Element( 93, 'Np', 'Neptunium'),   'Pu': Element( 94, 'Pu', 'Plutonium'),
    'Am': Element( 95, 'Am', 'Americium'),   'Cm': Element( 96, 'Cm', 'Curium'),
    'Bk': Element( 97, 'Bk', 'Berkelium'),   'Cf': Element( 98, 'Cf', 'Californium'),
    'Es': Element( 99, 'Es', 'Einsteinium'), 'Fm': Element(100, 'Fm', 'Fermium'),
    'Md': Element(101, 'Md', 'Mendelevium'), 'No': Element(102, 'No', 'Nobelium'),
    'Lr': Element(103, 'Lr', 'Lawrencium'),  'Rf': Element(104, 'Rf', 'Rutherfordium'),
    'Db': Element(105, 'Db', 'Dubnium'),     'Sg': Element(106, 'Sg', 'Seaborgium'),
    'Bh': Element(107, 'Bh', 'Bohrium'),     'Hs': Element(108, 'Hs', 'Hassium'),
    'Mt': Element(109, 'Mt', 'Meitnerium'),  'Ds': Element(110, 'Ds', 'Darmstadtium'),
    'Rg': Element(111, 'Rg', 'Roentgenium'), 'Cn': Element(112, 'Cn', 'Copernicium')
}