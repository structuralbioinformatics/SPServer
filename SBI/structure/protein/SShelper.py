import numpy as np
from collections        import Counter
from random             import randint
from ..atom             import AtomOfAminoAcid
from SBI.external.DSSP  import DSSPExe
from Arch               import Arch
from SecondaryStructure import SecondaryStructure


def calculate_dssp(pdb, tmppdb=None, tmpdssp=None, cleanfiles=True):

    prefix     = ".".join([pdb.globalID, str(randint(0, 100))])
    if tmppdb is None:
        tmppdb  = ".".join([prefix, 'pdb'])
    if tmpdssp is None:
        tmpdssp = ".".join([prefix, 'dssp'])

    newchain = pdb.duplicate(hetero=False, water=False)
    newchain.clean()
    fd = open(tmppdb, 'w')
    fd.write(newchain.PDB_format())
    fd.close()

    dsspexe = DSSPExe(tmppdb, tmpdssp, cleanfiles, cleanfiles)
    m = 0
    for n in range(len(pdb.aminoacids)):
        if pdb.aminoacids[n].single_letter == dsspexe.dsspdata[m].aminoacid:
            if pdb.aminoacids[n].has_full_backbone:
                pdb.aminoacids[n].dssp = dsspexe.dsspdata[m]
                step = 1 if m+1 < len(dsspexe.dsspdata) else 0
                m += step
            else:
                pdb.aminoacids[n].dssp = dsspexe.empty_dssp
        else:
            pdb.aminoacids[n].dssp = dsspexe.empty_dssp

    if dsspexe.gapped:
        template = list(pdb.gapped_protein_secondary_structure)
        for i in range(len(template)):
            if template[i] == 'x':
                pdb._gaps.append(i)


def locate_valid_secondarystructures(pdb):

    if not pdb.has_dssp:
        pdb.calculate_dssp()

    counter, start_post, reading_type, added = 0, None, 'C', False
    previous_aa = None

    for residue in pdb.aminoacids:
        ss, aa = residue.secondary_structure, residue.identifier
        if reading_type not in SecondaryStructure.structure_regions:
            counter = 0
            if ss in SecondaryStructure.structure_regions:
                start_post, reading_type, counter = aa, ss, 1
        else:
            if ss == reading_type:
                counter += 1
                if counter == SecondaryStructure.min_ss_length[ss]:
                    pdb._ss_valid.append(SecondaryStructure(ss, start_post))
                    added = True
            else:
                if added:
                    pdb._ss_valid[-1]._length = counter
                    pdb._ss_valid[-1]._endp   = previous_aa
                counter, start_post, reading_type, added = 0, aa, ss, False
                if ss in SecondaryStructure.structure_regions:
                    start_post, reading_type, counter = aa, ss, 1

        previous_aa = aa

    if len(pdb._ss_valid) > 0:
        if pdb._ss_valid[-1]._endp is None:
            pdb._ss_valid[-1]._length = counter
            pdb._ss_valid[-1]._endp   = previous_aa

    for ss in pdb._ss_valid:
        ss._struct = pdb.extract(ss._inip, ss._endp)
        ss.calculate_center_of_masses()


def calculate_archs(pdb, limit_internal_ss=100,
                    limit_distance=False, allowed_gaps=0):
    number = 1
    for i in range(len(pdb.secondary_structures)):
        for j in range(i+1, len(pdb.secondary_structures)):
            f1 = pdb._get_structure_array_coordinate(pdb.secondary_structures[i]._inip)
            f2 = pdb._get_structure_array_coordinate(pdb.secondary_structures[i]._endp)
            f3 = pdb._get_structure_array_coordinate(pdb.secondary_structures[j]._inip)
            f4 = pdb._get_structure_array_coordinate(pdb.secondary_structures[j]._endp)

            internal_ss_header = []
            if j-i-1 > limit_internal_ss:
                break
            else:
                ss_count = 2
                intss_types = ''
                for n in range(i+1, j):
                    internal_ss_header.append(pdb.secondary_structures[n].headerformat(ss_count, chr(65 + ss_count - 1)))
                    intss_types += pdb.secondary_structures[n]._sstype
                    ss_count += 1
            newss = pdb.extract(pdb.secondary_structures[i]._inip, pdb.secondary_structures[j]._endp)
            ss_template = Counter(list(newss.gapped_protein_secondary_structure))
            if 'x' in ss_template and ss_template['x'] > allowed_gaps:
                break
            gapcorrectini = 0
            gapcorrectend = 0
            if len(pdb._gaps) > 0:
                if f1 >= pdb._gaps[0]:
                    tmp_str       = pdb.extract(pdb.aminoacids[0].identifier, pdb.secondary_structures[i]._inip)
                    ss_template   = Counter(list(tmp_str.gapped_protein_secondary_structure))
                    gapcorrectini = ss_template['x']
                    if allowed_gaps == 0:
                        gapcorrectend = gapcorrectini
                    else:
                        tmp_str       = pdb.extract(pdb.aminoacids[0].identifier, pdb.secondary_structures[j]._endp)
                        ss_template   = Counter(list(tmp_str.gapped_protein_secondary_structure))
                        gapcorrectend = ss_template['x']

            secondstructurepair = Arch(pdb.globalID,
                                       pdb.secondary_structures[i],
                                       pdb.secondary_structures[j],
                                       j-i-1, f3 - f2 - 1, newss,
                                       pdb._torsionsCA[f1+gapcorrectini:f4+1+gapcorrectend],
                                       pdb._psiphi[(f1+gapcorrectini)*3:(f4+gapcorrectend)*3+6],
                                       number)
            # distance restriction only aplies to archs with internal structures
            if secondstructurepair.is_superarch and \
               (limit_distance and secondstructurepair.cartesian_distance > limit_distance):
                break

            secondstructurepair._inttxt = internal_ss_header
            secondstructurepair._inttyp = intss_types
            if secondstructurepair.is_superarch: pdb._superarchs.append(secondstructurepair)
            else:                                pdb._archs.append(secondstructurepair)
            number += 1


def calculate_torsionsCA(pdb, accuracy = 0):
    for i in range(0, len(pdb.aminoacids) - 3):
        aa0 = pdb.aminoacids[i]
        aa1 = pdb.aminoacids[i+1]
        aa2 = pdb.aminoacids[i+2]
        aa3 = pdb.aminoacids[i+3]

        ca0 = aa0.ca
        ca1 = aa1.ca
        ca2 = aa2.ca
        ca3 = aa3.ca

        if not isinstance(ca0, AtomOfAminoAcid) or not isinstance(ca1, AtomOfAminoAcid) or \
           not isinstance(ca2, AtomOfAminoAcid) or not isinstance(ca3, AtomOfAminoAcid):
            pdb._torsionsCA += '-'
            continue

        d01 = ca0.distance(ca1)
        v01 = np.divide(np.subtract(ca1.coordinates, ca0.coordinates), d01)
        d12 = ca1.distance(ca2)
        v12 = np.divide(np.subtract(ca2.coordinates, ca1.coordinates), d12)
        d23 = ca2.distance(ca3)
        v23 = np.divide(np.subtract(ca3.coordinates, ca2.coordinates), d23)

        c4  = -1 * np.sum(np.multiply(v01, v12))
        s4  = np.sqrt(1 - np.power(c4, 2))
        a4  = np.degrees(np.arctan(s4/c4))
        if a4 < 0: a4 += 180
        c5  = -1 * np.sum(np.multiply(v12, v23))
        s5  = np.sqrt(1 - np.power(c5, 2))
        a5  = np.degrees(np.arctan(s5/c5))
        if a5 < 0: a5 += 180

        u   = (v01[1]*v12[2])-(v01[2]*v12[1]), (v01[2]*v12[0])-(v01[0]*v12[2]), (v01[0]*v12[1])-(v01[1]*v12[0])
        v   = (v12[1]*v23[2])-(v12[2]*v23[1]), (v12[2]*v23[0])-(v12[0]*v23[2]), (v12[0]*v23[1])-(v12[1]*v23[0])
        c6  = np.sum(np.multiply(u, v))/(s4*s5)
        s6  = np.sum(np.multiply(v01, v))/(s4*s5)
        a6  = np.degrees(np.arctan(s6/c6))
        a7  = a6
        if a7 < 0 and s6 > 0: a6 += 180
        if a7 > 0 and c6 < 0: a6 -= 180
        if a6 < 0:            a6 += 360

        accuracy_desc = {'4':  ['A', 'B', 'C', 'D'],
                         '12': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']}
        if accuracy in accuracy_desc:
            number_of_parts = len(accuracy_desc[str(accuracy)])
            for j in range(number_of_parts):
                if a6 >= j*(360/number_of_parts) and a6 < (j+1)*(360/number_of_parts):
                    pdb._torsionsCA += accuracy_desc[str(accuracy)][j]
                    break
                if a6 == 360:
                    pdb._torsionsCA += accuracy_desc[str(accuracy)][-1]
        else:
            if a6 >= 0   and a6 < 90:   pdb._torsionsCA += 'A'
            if a6 >= 90  and a6 < 160:  pdb._torsionsCA += 'B'
            if a6 >= 160 and a6 < 260:  pdb._torsionsCA += 'C'
            if a6 >= 260 and a6 <= 360: pdb._torsionsCA += 'D'

    pdb._torsionsCA += '-'

    if len(pdb._gaps) > 0:
        template = list(pdb.gapped_protein_secondary_structure)
        tomodify = list(pdb._torsionsCA)
        for x in range(len(template)):
            if template[x] == 'x' and template[x-1] != 'x':
                tomodify[x-1] = '-'
                tomodify[x]   = '-'
            elif template[x-1] == 'x':
                if x < len(tomodify):
                    tomodify[x]   = '-'
                tomodify.insert(x+1, '-')
        pdb._torsionsCA = "".join(tomodify)


def calculate_psiphi(pdb):
    step  = 0
    phi   = []  # step 0 = C(n)   - N(n+1) - CA(n+1) - C(n+1)  : phi
    phi.append(0)
    psi   = 0   # step 1 = N(n-1) - CA(n)  - C(n)    - N(n)    : psi
    omega = []  # step 2 = CA(n)  - C(n)   - N(n+1)  - CA(n+1) : omega
    omega.append(0)
    m     = 0
    uncalculable = 0

    for i in range(3*(len(pdb.aminoacids)-1)):
        if step == 0:
            at0 = pdb.aminoacids[m].c
            at1 = pdb.aminoacids[m+1].n
            at2 = pdb.aminoacids[m+1].ca
            at3 = pdb.aminoacids[m+1].c
        elif step == 1:
            at0 = pdb.aminoacids[m].n
            at1 = pdb.aminoacids[m].ca
            at2 = pdb.aminoacids[m].c
            at3 = pdb.aminoacids[m+1].n
        else:
            at0 = pdb.aminoacids[m].ca
            at1 = pdb.aminoacids[m].c
            at2 = pdb.aminoacids[m+1].n
            at3 = pdb.aminoacids[m+1].ca

        if not isinstance(at0, AtomOfAminoAcid) or not isinstance(at1, AtomOfAminoAcid) or \
           not isinstance(at2, AtomOfAminoAcid) or not isinstance(at3, AtomOfAminoAcid) or \
           (step == 2 and uncalculable > 0):
            if step < 2:
                uncalculable += 1
            else:
                if uncalculable > 0:
                    pdb._psiphi += '---'
                    uncalculable = 0
            step += 1
            if step == 3:
                step = 0
                m   += 1
            continue

        d01 = at0.distance(at1)
        v01 = np.divide(np.subtract(at1.coordinates, at0.coordinates), d01)
        d12 = at1.distance(at2)
        v12 = np.divide(np.subtract(at2.coordinates, at1.coordinates), d12)
        d23 = at2.distance(at3)
        v23 = np.divide(np.subtract(at3.coordinates, at2.coordinates), d23)

        c4  = -1 * np.sum(np.multiply(v01, v12))
        s4  = np.sqrt(1 - np.power(c4, 2))
        c5  = -1 * np.sum(np.multiply(v12, v23))
        s5  = np.sqrt(1 - np.power(c5, 2))

        u   = (v01[1]*v12[2])-(v01[2]*v12[1]), (v01[2]*v12[0])-(v01[0]*v12[2]), (v01[0]*v12[1])-(v01[1]*v12[0])
        v   = (v12[1]*v23[2])-(v12[2]*v23[1]), (v12[2]*v23[0])-(v12[0]*v23[2]), (v12[0]*v23[1])-(v12[1]*v23[0])
        c6  = np.sum(np.multiply(u, v))/(s4*s5)
        s6  = np.sum(np.multiply(v01, v))/(s4*s5)
        a6  = np.degrees(np.arctan(s6/c6))
        a7  = a6
        if a7 < 0 and s6 > 0: a6 += 180
        if a7 > 0 and c6 < 0: a6 -= 180

        if step == 0: phi.append(a6)
        if step == 1: psi   = a6
        if step == 2:
            omega.append(a6)
            if omega[-2] >= -60 and omega[-2] <= 60 and len(omega) > 2:
                pdb._psiphi += 'CIS'  # cys-proline
            else:
                if phi[-2] >= -180 and phi[-2] < -20:
                    if phi[-2] >= -180 and phi[-2] < -140:
                        if   psi >= -180 and psi < -140: pdb._psiphi += 'E11'
                        elif psi >= -140 and psi < -100: pdb._psiphi += 'E12'
                        elif psi >= -100 and psi <  -60: pdb._psiphi += 'I13'
                        elif psi >=  -60 and psi <  -20: pdb._psiphi += 'N14'
                        elif psi >=  -20 and psi <   20: pdb._psiphi += 'N15'
                        elif psi >=   20 and psi <   60: pdb._psiphi += 'N16'
                        elif psi >=   60 and psi <  100: pdb._psiphi += 'E17'
                        elif psi >=  100 and psi <  140: pdb._psiphi += 'E18'
                        elif psi >=  140 and psi <= 180: pdb._psiphi += 'E19'
                    elif phi[-2] >= -140 and phi[-2] < -100:
                        if   psi >= -180 and psi < -140: pdb._psiphi += 'B21'
                        elif psi >= -140 and psi < -100: pdb._psiphi += 'F22'
                        elif psi >= -100 and psi <  -60: pdb._psiphi += 'I23'
                        elif psi >=  -60 and psi <  -20: pdb._psiphi += 'I24'
                        elif psi >=  -20 and psi <   20: pdb._psiphi += 'T25'
                        elif psi >=   20 and psi <   60: pdb._psiphi += 'T26'
                        elif psi >=   60 and psi <  100: pdb._psiphi += 'B27'
                        elif psi >=  100 and psi <  140: pdb._psiphi += 'B28'
                        elif psi >=  140 and psi <= 180: pdb._psiphi += 'B29'
                    elif phi[-2] >= -100 and phi[-2] < -60:
                        if   psi >= -180 and psi < -140: pdb._psiphi += 'B31'
                        elif psi >= -140 and psi < -100: pdb._psiphi += 'F32'
                        elif psi >= -100 and psi <  -60: pdb._psiphi += 'H33'
                        elif psi >=  -60 and psi <  -20: pdb._psiphi += 'H34'
                        elif psi >=  -20 and psi <   20: pdb._psiphi += 'T35'
                        elif psi >=   20 and psi <   60: pdb._psiphi += 'T36'
                        elif psi >=   60 and psi <  100: pdb._psiphi += 'B37'
                        elif psi >=  100 and psi <  140: pdb._psiphi += 'B38'
                        elif psi >=  140 and psi <= 180: pdb._psiphi += 'B39'
                    elif phi[-2] >= -60 and phi[-2] < -20:
                        if   psi >= -180 and psi < -140: pdb._psiphi += 'P41'
                        elif psi >= -140 and psi < -100: pdb._psiphi += 'F42'
                        elif psi >= -100 and psi <  -60: pdb._psiphi += 'H43'
                        elif psi >=  -60 and psi <  -20: pdb._psiphi += 'H44'
                        elif psi >=  -20 and psi <   20: pdb._psiphi += 'T45'
                        elif psi >=   20 and psi <   60: pdb._psiphi += 'T46'
                        elif psi >=   60 and psi <  100: pdb._psiphi += 'P47'
                        elif psi >=  100 and psi <  140: pdb._psiphi += 'P48'
                        elif psi >=  140 and psi <= 180: pdb._psiphi += 'P49'
                elif phi[-2] >= 20 and phi[-2] <= 180:
                    if phi[-2] >= 20 and phi[-2] < 60:
                        if   psi >= -180 and psi < -140: pdb._psiphi += 'G61'
                        elif psi >= -140 and psi < -100: pdb._psiphi += 'G62'
                        elif psi >= -100 and psi <  -60: pdb._psiphi += 'G63'
                        elif psi >=  -60 and psi <  -20: pdb._psiphi += 'U64'
                        elif psi >=  -20 and psi <   20: pdb._psiphi += 'U65'
                        elif psi >=   20 and psi <   60: pdb._psiphi += 'L66'
                        elif psi >=   60 and psi <  100: pdb._psiphi += 'L67'
                        elif psi >=  100 and psi <  140: pdb._psiphi += 'M68'
                        elif psi >=  140 and psi <= 180: pdb._psiphi += 'M69'
                    elif phi[-2] >= 60 and phi[-2] < 100:
                        if   psi >= -180 and psi < -140: pdb._psiphi += 'G71'
                        elif psi >= -140 and psi < -100: pdb._psiphi += 'G72'
                        elif psi >= -100 and psi <  -60: pdb._psiphi += 'G73'
                        elif psi >=  -60 and psi <  -20: pdb._psiphi += 'U74'
                        elif psi >=  -20 and psi <   20: pdb._psiphi += 'U75'
                        elif psi >=   20 and psi <   60: pdb._psiphi += 'L76'
                        elif psi >=   60 and psi <  100: pdb._psiphi += 'L77'
                        elif psi >=  100 and psi <  140: pdb._psiphi += 'M78'
                        elif psi >=  140 and psi <= 180: pdb._psiphi += 'M79'
                    elif phi[-2] >= 100 and phi[-2] < 140:
                        if   psi >= -180 and psi < -140: pdb._psiphi += 'G81'
                        elif psi >= -140 and psi < -100: pdb._psiphi += 'G82'
                        elif psi >= -100 and psi <  -60: pdb._psiphi += 'G83'
                        elif psi >= -60  and psi <  -20: pdb._psiphi += 'U84'
                        elif psi >= -20  and psi <   20: pdb._psiphi += 'U85'
                        elif psi >= 20   and psi <   60: pdb._psiphi += 'S86'
                        elif psi >= 60   and psi <  100: pdb._psiphi += 'S87'
                        elif psi >= 100  and psi <  140: pdb._psiphi += 'S88'
                        elif psi >=140   and psi <= 180: pdb._psiphi += 'S89'
                    elif phi[-2] >= 140 and phi[-2] <= 180:
                        if   psi >= -180 and psi < -140: pdb._psiphi += 'E91'
                        elif psi >= -140 and psi < -100: pdb._psiphi += 'E92'
                        elif psi >= -100 and psi <  -60: pdb._psiphi += 'I93'
                        elif psi >=  -60 and psi <  -20: pdb._psiphi += 'N94'
                        elif psi >=  -20 and psi <   20: pdb._psiphi += 'N95'
                        elif psi >=   20 and psi <   60: pdb._psiphi += 'N96'
                        elif psi >=   60 and psi <  100: pdb._psiphi += 'E97'
                        elif psi >=  100 and psi <  140: pdb._psiphi += 'E98'
                        elif psi >=  140 and psi <= 180: pdb._psiphi += 'E99'
                elif phi[-2] >= -20 and phi[-2] <20:
                    if   psi >= -180 and psi < -140: pdb._psiphi += 'O51'
                    elif psi >= -140 and psi < -100: pdb._psiphi += 'O52'
                    elif psi >= -100 and psi <  -60: pdb._psiphi += '*53'
                    elif psi >=  -60 and psi <  -20: pdb._psiphi += '*54'
                    elif psi >=  -20 and psi <   20: pdb._psiphi += '*55'
                    elif psi >=   20 and psi <   60: pdb._psiphi += '*56'
                    elif psi >=   60 and psi <  100: pdb._psiphi += '*57'
                    elif psi >=  100 and psi <  140: pdb._psiphi += 'O58'
                    elif psi >=  140 and psi <= 180: pdb._psiphi += 'O59'
        step += 1
        if step == 3:
            step = 0
            m   += 1

    pdb._psiphi += '---'

    if len(pdb._gaps) > 0:
        template = list(pdb.gapped_protein_secondary_structure)
        tomodify = [pdb._psiphi[i:i+3] for i in range(0, len(pdb._psiphi), 3)]
        for x in range(len(template)):
            if template[x] == 'x' and template[x-1] != 'x':
                tomodify[x-1] = '---'
                tomodify.insert(x, '---')
            elif template[x] == 'x':
                tomodify.insert(x, '---')
        pdb._psiphi = "".join(tomodify)
