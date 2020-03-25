from SBI.beans.IndexedNum import IndexedNum
from SBI.data import aminoacids3to1, nucleic2to1
import time
import re


def process_HEADER_line(line):
    '''
    COLUMNS  DATA TYPE    FIELD           DEFINITION
    ---------------------------------------------------------------------------
     1 -  6  Record name  "HEADER"
    11 - 50  String(40)   classification  Classifies the molecule(s)
    51 - 59  Date         depDate         Deposition date.
                                          his is the date the coordinates were
                                          received by the PDB
    63 - 66  IDcode       idCode          This identifier is unique in the PDB
    '''
    data = []
    data.append(line[10:50].strip())
    inf, outf = "%d-%b-%y", "%Y-%m-%d"
    if line[50:60].strip() != '':
        tme = time.strftime(outf, time.strptime(line[50:61].strip(), inf))
    else:
        tme = ''
    data.append(tme)
    data.append(line[62:67].strip())
    return data

def process_OBSLTE_line(line):
    '''
    COLUMNS     DATA TYPE       FIELD          DEFINITION
    -------------------------------------------------------------------------------
     1 -  6     Record name     "OBSLTE"

     9 - 10     Continuation    continuation   Allows concatenation of multiple
                                               records.

    12 - 20     Date            repDate        Date that this entry was replaced.

    22 - 25     IDcode          idCode         ID code of this entry.

    32 - 35     IDcode          rIdCode        ID code of entry that replaced
                                               this one.

    37 - 40     IDcode          rIdCode        ID code of entry that replaced
                                               this one.

    42 - 45     IDcode          rIdCode        ID code of entry that replaced
                                               this one.

    47 - 50     IDcode          rIdCode        ID code of entry that replaced
                                               this one.

    52 - 55     IDcode          rIdCode        ID code of entry that replaced
                                               this one.

    57 - 60     IDcode          rIdCode        ID code of entry that replaced
                                               this one.

    62 - 65     IDcode          rIdCode        ID code of entry that replaced
                                               this one.

    67 - 70     IDcode          rIdCode        ID code of entry that replaced
                                               this one.
    '''
    data = []
    inf, outf = "%d-%b-%y", "%Y-%m-%d"
    if line[11:21].strip() != '':
        tme = time.strftime(outf, time.strptime(line[11:21].strip(), inf))
    else:
        tme = ''
    data.append(tme)
    data.extend(line[31:].strip().split())
    return data

def process_TITLE_line(line):
    '''
    COLUMNS  DATA TYPE     FIELD         DEFINITION
    --------------------------------------------------------------------------
     1 -  6  Record name   "TITLE "
     9 - 10  Continuation  continuation  Concatenation of multiple records.
    11 - 70  String        title         Title of the experiment.
    '''
    return line[10:].strip()


def process_EXPERIMENT_line(line):
    '''
    COLUMNS       DATA TYPE      FIELD         DEFINITION
    -------------------------------------------------------
     1 -  6       Record name    "EXPDTA"
     9 - 10       Continuation   continuation  Allows concatenation
                                               of multiple records
    11 - 70       SList          technique     The experimental technique(s)
                                               with optional comment describing
                                               the sample or experiment.
    '''
    return line[10:].split(';')[0].strip()


def process_RESOLUTION_line(line):
    '''
    COLUMNS    DATA TYPE      FIELD              DEFINITION
    -------------------------------------------------------
     1 - 6     Record name    "REMARK"
    10         LString(1)     "2"
    12 - 22    LString(11)    "RESOLUTION."
    23 - 27    Real(5.2)      resolution         Resolution.
    29 - 38    LString(10)    "ANGSTROMS."
    '''
    data = line[23:30].strip().split()[0]
    return None if data == 'NULL' or data == 'NOT' else data


def process_RFACTOR_line(line):
    if line.strip().endswith('NULL'):
        return None
    return line.split(':')[-1].strip()


def process_FREER_line(line):
    if line.strip().endswith('NULL'):
        return None
    return line.split(':')[-1].strip()


def process_SUPERSEEDED_line(line):
    '''
    COLUMNS        DATA TYPE     FIELD         DEFINITION
    -------------------------------------------------------------------------
     1 -  6        Record name   "SPRSDE"
     9 - 10        Continuation  continuation  Allows for multiple ID codes.
    12 - 20        Date          sprsdeDate    Date it superseded others
    22 - 25        IDcode        idCode        ID code of this entry.
    32 - 35        IDcode        sIdCode       ID code of a superseded entry.
    37 - 40        IDcode        sIdCode       ID code of a superseded entry.
    42 - 45        IDcode        sIdCode       ID code of a superseded entry.
    47 - 50        IDcode        sIdCode       ID code of a superseded entry.
    52 - 55        IDcode        sIdCode       ID code of a superseded entry.
    57 - 60        IDcode        sIdCode       ID code of a superseded entry.
    62 - 65        IDcode        sIdCode       ID code of a superseded entry.
    67 - 70        IDcode        sIdCode       ID code of a superseded entry.
    72 - 75        IDcode        sIdCode       ID code of a superseded entry.
    '''
    return line[10:].strip().split()[2:]


def process_MOLKEY(line):
    key = line.split(':')[-1].strip()
    return re.sub(';', '', key)


def process_COMPND_line(line):
    '''
    COLUMNS  DATA TYPE       FIELD         DEFINITION
    --------------------------------------------------------------------------
     1 -  6  Record name     "COMPND"
     8 - 10  Continuation    continuation  Concatenation of multiple records.
    11 - 80  Specification   compound      Description of molecular components.
             list
    '''
    return line[10:].strip()


def process_SOURCE_line(line):
    '''
    COLUMNS DATA  TYPE     FIELD          DEFINITION
    -------------------------------------------------------------------------
     1 -  6 Record name    "SOURCE"
     8 - 10 Continuation   continuation   Concatenation of multiple records.
    11 - 79 Specification  srcName        Identifies the source of the
            List                          macromolecule in a token.
    '''
    return line[10:].strip()


def process_KEYWRD_line(line):
    '''
        COLUMNS       DATA  TYPE     FIELD         DEFINITION
    -------------------------------------------------------------------
     1 -  6       Record name    "KEYWDS"
     9 - 10       Continuation   continuation  Concatenation of records.
    11 - 79       List           keywds        Comma-separated list
    '''
    return line[10:].rstrip()


def process_SITE_IDENTIFIER_line(line):
    '''
             1         2         3         4         5         6
    123456789012345678901234567890123456789012345678901234567890
    REMARK 800
    REMARK 800 SITE
    REMARK 800 SITE_IDENTIFIER: FREE TEXT GOES HERE.
    REMARK 800 EVIDENCE_CODE: (AUTHOR or SOFTWARE or UNKNOWN)
    REMARK 800 SITE_DESCRIPTION: FREE TEXT GOES HERE.
    '''
    return line[10:].split(':')[1].strip()


def process_REMARK800_line(line):
    return line[10:].strip()


def process_SITE_lines(line):
    '''
             1         2         3         4         5         6
    123456789012345678901234567890123456789012345678901234567890
    SITE     1 AC1  3 HIS A  94  HIS A  96  HIS A 119
    SITE     1 AC2  5 ASN A  62  GLY A  63  HIS A  64  HOH A 328
    SITE     2 AC2  5 HOH A 634
    SITE     1 AC3  5 GLN A 136  GLN A 137  PRO A 138  GLU A 205
    SITE     2 AC3  5 CYS A 206
    SITE     1 AC4 11 HIS A  64  HIS A  94  HIS A  96  HIS A 119
    SITE     2 AC4 11 LEU A 198  THR A 199  THR A 200  TRP A 209
    SITE     3 AC4 11 HOH A 572  HOH A 582  HOH A 635
    '''
    data   = [line[10:14].strip()]
    line   = line[18:].rstrip()
    i = 0
    while i + 9 <= len(line):
        resi = line[i:i + 11]
        rest = resi[:3].strip()
        resp = resi[3:].strip()
        data.append(tuple([rest, resp]))
        i += 11
    return data


def process_HETEROATOM_lines(line):
    '''
    There are 3 lines that define an heteroatom:

    COLUMNS     DATA TYPE     FIELD         DEFINITION
    ------------------------------------------------------
     1 -  6     Record name   "HET      "
     8 - 10     LString(3)    hetID         Het identifier, right-justified.
    13          Character     ChainID       Chain identifier.
    14 - 17     Integer       seqNum        Sequence number.
    18          AChar         iCode         Insertion code.
    21 - 25     Integer       numHetAtoms   Number of HETATM records for the
                                            group present in the entry.
    31 - 70     String        text          Text describing Het group.

             1         2         3         4         5         6         7
    1234567890123456789012345678901234567890123456789012345678901234567890
    HET    TRS    975       8
    HET    STA  I   4      25     PART_OF: HIV INHIBITOR;

    HET    FUC  Y   1      10     PART_OF: NONOATE COMPLEX; L-FUCOSE
    HET    GAL  Y   2      11     PART_OF: NONOATE COMPLEX
    HET    NAG  Y   3      15
    HET    FUC  Y   4      10
    HET    NON  Y   5      12

    HET    UNX  A 161       1     PSEUDO CARBON ATOM OF UNKNOWN LIGAND
    HET    UNX  A 162       1     PSEUDO CARBON ATOM OF UNKNOWN LIGAND
    HET    UNX  A 163       1     PSEUDO CARBON ATOM OF UNKNOWN LIGAND

    COLUMNS     DATA TYPE       FIELD          DEFINITION
    -------------------------------------------------------
     1 -  6     Record name     "HETNAM"
     9 - 10     Continuation    continuation   Allows concatenation of
                                               multiple records.
    12 - 14     LString(3)      hetID          Het identifier,
                                               right-justified.
    16 - 70     String          text           Chemical name.

             1         2         3         4         5         6         7
    1234567890123456789012345678901234567890123456789012345678901234567890
    HETNAM     GLC GLUCOSE
    HETNAM     SAD BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE
    HETNAM  2  SAD DINUCLEOTIDE

    HETNAM     UNX UNKNOWN ATOM OR ION
    HETNAM     UNL UNKNOWN LIGAND

    HETNAM     CYE 45-(3-AMINOPROPYL)-5,11,22,28,34-PENTAMETHYL-3,9,15,
    HETNAM   2 CYE  20,26,32,38,43-OCTAOXO-2,5,8,14,19,22,25,28,31,34,37,
    HETNAM   3 CYE  42,45,48-TETRADECAAZA-11-AZONIAHEPTACYCLO[42.2.1.1~4,
    HETNAM   4 CYE  7~.1~10,13~.1~21,24~.1~27,30~.1~33,36~]DOPENTACONTA-
    HETNAM   5 CYE  1(46),4(52),6,10(51),12,21(50),23,27(49),29,33(48),35,
    HETNAM   6 CYE  44(47)-DODECAENE

    COLUMNS     DATA TYPE       FIELD          DEFINITION
    -------------------------------------------------------
     1 -  6     Record name     "FORMUL"
     9 - 10     Integer         compNum        Component number.
    13 - 15     LString(3)      hetID          Het identifier.
    17 - 18     Integer         continuation   Continuation number.
    19          Character       asterisk       "*" for water.
    20 - 70     String          text           Chemical formula.

             1         2         3         4         5         6         7
    1234567890123456789012345678901234567890123456789012345678901234567890
    FORMUL   2  SO4    2(O4 S1 2-)
    FORMUL   3  GLC    C6 H12 O6
    FORMUL   3  FOL    2(C19 H17 N7 O6 2-)
    FORMUL   4   CL    2(CL1 1-)
    FORMUL   5   CA    CA1 2+

    FORMUL   1  ACE    C2 H3 O1
    FORMUL   2  ACE    C2 H3 O1

    FORMUL   8  HOH   *463(H2 O1)
    '''
    data = []
    if line.startswith("HET   "):
        data.append(line[5:].strip().split()[0])
        data.append(line[11:17].strip())
    elif line.startswith("HETNAM"):
        data.append(line[10:14].strip())
        data.append(line[15:].strip())
    elif line.startswith("FORMUL"):
        data.append(line[11:15].strip())
        formula = re.sub(r'^\d+\(', '', line[16:].strip())
        formula = re.sub(r'\)\Z', '', formula)
        data.append(formula)
    return data


def process_DBREF_line(line):
    '''
    The line format is such as:

    COLUMNS    DATA TYPE    FIELD          DEFINITION
    ----------------------------------------------------------------
     1 - 6     Record name  "DBREF "
     8 - 11    IDcode       idCode         ID code of this entry.
    13         Character    chainID        Chain identifier.
    15 - 18    Integer      seqBegin       Initial sequence number
                                           of the PDB sequence segment.
    19         AChar        insertBegin    Initial insertion code
                                           of the PDB sequence segment.
    21 - 24   Integer       seqEnd         Ending sequence number
                                           of the PDB sequence segment.
    25         AChar        insertEnd      Ending insertion code
                                           of the PDB sequence segment.
    27 - 32    LString      database       Sequence database name.
    34 - 41    LString      dbAccession    Sequence database accession code.
    43 - 54    LString      dbIdCode       Sequence database
                                           identification code.
    56 - 60    Integer      dbseqBegin     Initial sequence number of the
                                           database seqment.
    61         AChar        idbnsBeg       Insertion code of initial residue
                                           of the segment, if PDB is the
                                           reference.
    63 - 67    Integer      dbseqEnd       Ending sequence number of the
                                           database segment.
    68         AChar        dbinsEnd       Insertion code of the ending
                                           residue of the segment, if PDB is
                                           the reference.

    For example:

             1         2         3         4         5         6         7
    1234567890123456789012345678901234567890123456789012345678901234567890
    DBREF  1ABC B    1B   36  PDB    1ABC     1ABC             1B    36

    DBREF  3AKY      3   220  SWS    P07170   KAD1_YEAST       5    222

    DBREF  1HAN      2   288  GB     397884   X66122           1    287

    DBREF  3HSV A    1    92  SWS    P22121   HSF_KLULA      193    284
    DBREF  3HSV B    1    92  SWS    P22121   HSF_KLULA      193    284

    DBREF  1ARL      1   307  SWS    P00730   CBPA_BOVIN     111    417

    DBREF  249D A    1    12  NDB    BDL070   BDL070           1     12
    DBREF  249D B   13    24  NDB    BDL070   BDL070          13     24
    DBREF  249D C   26    36  NDB    BDL070   BDL070          26     36
    DBREF  249D D   37    48  NDB    BDL070   BDL070          37     48

    Detailed info @ http://www.wwpdb.org/documentation/format23/sect3.html
    '''
    data = []
    data.append(line[7:11].strip())
    data.append(fix_chain(line[12:13].strip()))
    ini = line[14:19].strip()
    end = line[20:25].strip()
    data.append(IndexedNum(ini) if ini != '' else '')
    data.append(IndexedNum(end) if end != '' else '')
    data.append(line[25:29].strip())
    data.append(line[30:41].strip())
    return data


def process_SECONDARYSTRUCTURE_lines(line):
    '''
    COLUMNS  DATA  TYPE     FIELD         DEFINITION
    ---------------------------------------------------------------------------
     1 -  6  Record name    "HELIX "
     8 - 10  Integer        serNum        Serial number of the helix. Starts
                                          at 1  and increases incrementally.
    12 - 14  LString(3)     helixID       Helix  identifier.
    16 - 18  Residue name   initResName   Name of the initial residue.
    20       Character      initChainID   Chain identifier
    22 - 25  Integer        initSeqNum    Number of the initial residue.
    26       AChar          initICode     Code of the initial residue.
    28 - 30  Residue  name  endResName    Name of the terminal residue
    32       Character      endChainID    Chain identifier
    34 - 37  Integer        endSeqNum     Number of the terminal residue.
    38       AChar          endICode      Code of the terminal residue.
    39 - 40  Integer        helixClass    Helix class (see below).
    41 - 70  String         comment       Comment about this helix.
    72 - 76  Integer        length        Length of this helix.

    COLUMNS  DATA  TYPE     FIELD          DEFINITION
    --------------------------------------------------------------------------
     1 -  6  Record name   "SHEET "
     8 - 10  Integer       strand         Strand number
    12 - 14  LString(3)    sheetID        Sheet  identifier.
    15 - 16  Integer       numStrands     Number  of strands in sheet.
    18 - 20  Residue name  initResName    Residue  name of initial residue.
    22       Character     initChainID    Chain identifier of initial residue
                                          in strand.
    23 - 26  Integer       initSeqNum     Sequence number of initial residue
                                          in strand.
    27       AChar         initICode      Insertion code of initial residue
                                          in  strand.
    29 - 31  Residue name  endResName     Residue name of terminal residue.
    33       Character     endChainID     Chain identifier of terminal residue.
    34 - 37  Integer       endSeqNum      Sequence number of terminal residue.
    38       AChar         endICode       Insertion code of terminal residue.
    39 - 40  Integer       sense          Sense of strand.
                                          0 if first strand, 1 if  parallel
                                          and -1 if anti-parallel.
    42 - 45  Atom          curAtom        Atom name in current strand.
    46 - 48  Residue name  curResName     Residue name in current strand
    50       Character     curChainId     Chain identifier in
                                          current strand.
    51 - 54  Integer       curResSeq      Residue sequence number
                                          in current strand.
    55       AChar         curICode       Insertion code in
                                          current strand.
    57 - 60  Atom          prevAtom       Atom name in previous strand.
    61 - 63  Residue name  prevResName    Registration.  Residue name in
                                          previous strand.
    65       Character     prevChainId    Registration.  Chain identifier in
                                          previous  strand.
    66 - 69  Integer       prevResSeq     Registration. Residue sequence number
                                          in previous strand.
    70       AChar         prevICode      Registration.  Insertion code in
                                          previous strand.
    '''
    return line.strip()


def process_MATRIX_lines(line):
    '''
    COLUMNS        DATA  TYPE    FIELD         DEFINITION
    --------------------------------------------------------------------------
     1 -  6        Record name   "MTRIXn"      n=1, 2, or 3
     8 - 10        Integer       serial        Serial number.
    11 - 20        Real(10.6)    m[n][1]       Mn1
    21 - 30        Real(10.6)    m[n][2]       Mn2
    31 - 40        Real(10.6)    m[n][3]       Mn3
    46 - 55        Real(10.5)    v[n]          Vn
    60             Integer       iGiven
    '''
    # Matrix data can and MUST be gathered not by position
    # but by spliting the line
    # Why? 1M4X has more than 1000 matrices and shifts the position
    #      1QKP is wrongly formated and shifts the position
    # But... then there is 1RHV... obsolete, yes but...
    #      Check the number of segments to fix.
    l = line.replace('REMARK 290', '').replace('REMARK 350', '').strip()
    matrix_value    = re.split('\s+', l)
    if len(matrix_value) < 6:
        matrix_value = [l[:6].strip(),    l[7:10].strip(),
                        l[10:20].strip(), l[20:30].strip(),
                        l[30:40].strip(), l[45:55].strip()]
    matrix_value[0] = matrix_value[0][-1]
    matrix_value.pop(1)
    return matrix_value

def process_REMARK350_lines(line):
    return  line.strip().split(':')[-1]

def process_SEQRES_lines(line):
    '''
    COLUMNS        DATA TYPE       FIELD         DEFINITION
    ---------------------------------------------------------------------------
     1 -  6  Record name     "SEQRES"
     9 - 10  Integer         serNum        Serial number of the SEQRES record
                                           for the current chain.  Starts at 1
                                           and increments by one each line.
                                           Reset to 1 for each chain.
    12       Character       chainID       Chain identifier.  This may be any
                                           single legal character, including a
                                           blank which is used if there is
                                           only one chain.
    14 - 17  Integer         numRes        Number of residues in the chain.
                                           This value is repeated on every
                                           record.
    20 - 22  Residue name    resName       Residue name.
    24 - 26  Residue name    resName       Residue name.
    28 - 30  Residue name    resName       Residue name.
    32 - 34  Residue name    resName       Residue name.
    36 - 38  Residue name    resName       Residue name.
    40 - 42  Residue name    resName       Residue name.
    44 - 46  Residue name    resName       Residue name.
    48 - 50  Residue name    resName       Residue name.
    52 - 54  Residue name    resName       Residue name.
    56 - 58  Residue name    resName       Residue name.
    60 - 62  Residue name    resName       Residue name.
    64 - 66  Residue name    resName       Residue name.
    68 - 70  Residue name    resName       Residue name
    '''
    chain = fix_chain(line[11])
    data  = line[19:].strip().split()
    for i in range(len(data)):
        pos = data[i]
        if pos in aminoacids3to1:
            data[i] = aminoacids3to1[pos]
        elif pos in nucleic2to1:
            data[i] = nucleic2to1[pos]
        else:
            data[i] = 'X'
    return [chain, ''.join(data)]

def fix_chain(chain):
    '''
    Sometimes old PDBs have no chain, transform to chain A
    '''
    if chain.strip() == '' or chain.strip().upper() == 'NULL':
        return 'A'
    else:
        return chain.strip()
