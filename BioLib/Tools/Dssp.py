import os, sys, subprocess

def get_dssp_exec():
    '''
    Gets the DSSP executable
    '''
    try:
        dssp_path = os.environ['DSSP_PATH']
    except:
        sys.stderr.write('The $DSSP_PATH environment variable is not defined, trying local library dssp\n')
        return None
    return os.path.join(dssp_path, 'dsspcmbi')

def execute_dssp(structure):
    '''
    Executes dssp and sets secondary structure (ss) and accesible surface area (acc) to each residue of the given structure
    '''
    dssp_exec = get_dssp_exec()
    if not dssp_exec:
        dssp_exec = os.path.abspath(os.path.join(os.path.dirname(__file__), 'dssp'))
    if not os.path.exists(dssp_exec):
        sys.stderr.write('Cannot call dssp, skipping...\n')
        return None

    p = subprocess.Popen([dssp_exec, '--'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    dssp_out, dssp_err = p.communicate(str(structure))

    read = False
    for line in dssp_out.split('\n'):
        if line.startswith('  #  RESIDUE AA STRUCTURE BP1 BP2  ACC'):
            read = True
            continue
        if read and line != '' and line[13] != '!':
            residue_num = int(line[6:10].strip())
            ss = line[16:17]
            acc = int(line[35:38].strip())
            if ss != "H" and ss != "E":
                ss = "C"
            residue = structure.get_residue_by_num(residue_num)
            residue.set_ss(ss)
            residue.set_acc(acc)