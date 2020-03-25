import BioLib.Structure.PDB
import BioLib.Docking.FTDock

linkerPotential = {'2': (4.274, 11.645), 
                   '3': (5.558, 13.938),
                   '4': (6.651, 15.372),
                   '5': (7.088, 18.662), 
                   '6': (6.631, 22.224),
                   '7': (7.517, 24.180), 
                   '8': (8.312, 25.594), 
                   '9': (9.723, 29.490),
                   '10': (9.666, 33.388),
                   '11': (9.677, 36.293), 
                   '12': (8.461, 38.853), 
                   '13': (10.385, 38.106),
                   '14': (7.937, 40.773),
                   '15': (5.955, 53.153),
                  }

def decoy_linker_energy(static, mobile, ftdock, decoy, nativeMobile = None):
    '''
    Returns 1 if the decoy is compatible with linker statistics, it also returns
    the rmsd with the native solution if nativeMobile is provided
    '''
    staticDecoy = ftdock.get_static_structure(static)
    mobileDecoy = ftdock.get_mobile_structure(static, mobile, decoy)
    linkerDistance = staticDecoy.get_residues()[-1].get_ca_distance(mobileDecoy.get_residues()[0])
    linkerLength = mobileDecoy.get_residues()[0].get_num() - staticDecoy.get_residues()[-1].get_num() -1
    if linkerLength > 15:
         linkerLength = 15
    if linkerPotential[str(linkerLength)][1] > linkerDistance > linkerPotential[str(linkerLength)][0]:
         compatible = 1
    else:
         compatible = 0
    if nativeMobile:
        rmsd = ftdock.get_ligand_RMSD(static, mobile, nativeMobile, decoy)
        print str(decoy.get_id())+'\t'+str(rmsd)+'\t'+str(compatible)
        return compatible, rmsd
    return compatible

def rank_linker_enegy(static, mobile, ftdock, nativeMobile = None):
    '''
    Returns the compatibility of each decoy in a ftdock experiment
    '''
    linkerEnergyDecoy = {}
    for decoy in ftdock.get_decoys():
        linkerEnergyDecoy[decoy.get_id()] = decoy_linker_energy(static, mobile, ftdock, decoy, nativeMobile)
    return linkerEnergyDecoy
