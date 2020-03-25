def _sequencer(pdb, mode):
    seq = ''
    idx = []
    for aa in range(len(pdb.aminoacids)):
        if aa == 0:
            if mode == 'seq':
                seq += pdb.aminoacids[aa].single_letter
            if mode == 'str':
                seq += pdb.aminoacids[aa].secondary_structure
            if mode == 'exp':
                seq += pdb.aminoacids[aa].accessibilitycoded
            if mode == 'idx':
                idx.append(pdb.aminoacids[aa].identifier)
        else:
            if pdb.aminoacids[aa].follows(pdb.aminoacids[aa - 1]):
                if mode == 'seq':
                    seq += pdb.aminoacids[aa].single_letter
                if mode == 'str':
                    seq += pdb.aminoacids[aa].secondary_structure
                if mode == 'exp':
                    seq += pdb.aminoacids[aa].accessibilitycoded
                if mode == 'idx':
                    idx.append(pdb.aminoacids[aa].identifier)
            else:
                id_distance = -1 * \
                    pdb.aminoacids[aa].identifier_distance(
                        pdb.aminoacids[aa - 1])
                if id_distance > 1:
                    for x in range(id_distance - 1):
                        if mode == 'idx':
                            idx.append('X')
                        else:
                            seq += 'x'
                else:
                    if mode == 'idx':
                        idx.append('X')
                    else:
                        seq += 'x'
                if mode == 'seq':
                    seq += pdb.aminoacids[aa].single_letter
                if mode == 'str':
                    seq += pdb.aminoacids[aa].secondary_structure
                if mode == 'exp':
                    seq += pdb.aminoacids[aa].accessibilitycoded
                if mode == 'idx':
                    idx.append(pdb.aminoacids[aa].identifier)

    if mode == 'idx':
        return ";".join(idx)
    else:
        return seq
