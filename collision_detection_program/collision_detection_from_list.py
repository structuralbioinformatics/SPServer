#!/usr/bin/python

import sys,re
import numpy as np
from numpy import linalg as la 
import os
import math
import shlex,subprocess
import codecs

# Personal libraries
import SBI.structure.PDB as PDB

def interface(complexAB,moleculeAB):
    data = []
    com = np.array(complexAB)
    mol = np.array(moleculeAB)
    row = np.array(complexAB).shape[0]
    i = 0
    j = 0
    condition_i = 1 
    condition_j = 1

    for eachrow in range(1,row):
        if not (float (mol[eachrow,1]) == float (com[eachrow,1])):
            if (float (com[eachrow,1])):
                data.append(mol[eachrow,0:6])
                if ((int (mol[eachrow,0]) > i) & (condition_i)):
                    i = int (mol[eachrow,0])
                else : condition_i = 0 

        if ((int(mol[eachrow,0]) > j) & (condition_j)):
            j = int (mol[eachrow,0])
        else : condition_j = 0

    return np.array(data),i,j   
    
def normal_calculation(data,i,j,normal_vector):
    normal = []
    center_A = []
    center_B = []
    rectify = []
    resultant_molA = []
    resultant_molB = []
    mol = normal_vector
    for eachrow in range(1,len(mol)):
        if (int (mol[eachrow][1]) == 1):
            rectify.append(mol[eachrow])
        else : continue
    
    length = len(rectify)
    end = len(data)
    start_with_zero = 0
    Switch = 1 
    Switch_A = 1

    for row in range(0,end):
        for lne in range(start_with_zero,length):

            if (int (data[row,0]) == int (rectify[lne][0])):
                Switch = 0
                normal = rectify[lne][-7:]
                normal.append(rectify[lne][0])

                if (Switch_A == 1):
                    resultant_molA.append(normal)
                    center_A.append(data[row])
                else :
                    resultant_molB.append(normal)
                    center_B.append(data[row])

            elif (Switch == 0):
                if (int (data[row,0]) == i):
                    Switch_A = 0
                Switch = 1              
                start_with_zero = lne
                break

            else : continue 

    return np.array(resultant_molA),np.array(resultant_molB),np.array(center_A),np.array(center_B),j



def corresponding_molecule(normal_A,normal_B,center_A,center_B,k,pdb,Chain):
    atoms = []
    penetrating_res = []
    position_vec_AB = []
    atom_A = np.array(normal_A[:,7],dtype = int)
    atom_B = np.array(normal_B[:,7],dtype = int)
    pdb_AB = pdb
    example = np.matrix(normal_B[:,0])
    #A = np.array(center_A[:,2:5],np.float)
    #B = np.array(center_B[:,2:5],np.float)
    A = np.array(normal_A[:,0:3],np.float)
    B = np.array(normal_B[:,0:3],np.float)
    magnitude_A = (np.apply_along_axis(np.linalg.norm, 1, np.array(normal_A[:,4:7],np.float)))
    pdb = np.array(pdb_AB,dtype = 'string')
    temp = atom_A[0]
    uniq = []
    uniq = np.array(uniq,dtype = int)
    Area = 0
    condition = 1
    ligand_atoms = []
    ligand_atoms = np.array(ligand_atoms,dtype = int)
    Triangle = []

    for row in range(0,len(normal_A)):
        magnitude_distance = (np.apply_along_axis(np.linalg.norm,1,(np.cross(np.matrix(np.subtract(B,A[row,0:3])),np.matrix(normal_A[row,4:7],np.float)))))
        distance = np.array(np.divide(magnitude_distance,magnitude_A[row]))
        values = np.array(np.where(distance < 0.3))
        if not values.size: 
            continue
        Each_A = np.matrix(A[row,0:3],np.float)
        Each_B = np.matrix(B[values,0:3],np.float)
        position_vec_AB = np.subtract(Each_A,Each_B)
        normal_AB = np.apply_along_axis(np.linalg.norm,1,position_vec_AB)
        posit_unit_vec_AB = np.divide(position_vec_AB,normal_AB.reshape(len(normal_AB),1))
        normalA = np.matrix(normal_A[row,4:7],np.float)
        normalB = np.matrix(normal_B[values,4:7],np.float)
        dot_product_AB = np.append(np.sum(np.multiply(normalA,-position_vec_AB),axis=1),np.sum(np.multiply(normalB,position_vec_AB),axis=1),axis=1) # Numpy array of dot products
        s = np.array(np.where(dot_product_AB < 0))[0] # Get rows that contain negative values
        #rep_el = np.array(s[np.diff(s) == 0])
        rep_el = s[(np.where(np.diff(s) < 1))[0]]
        column = values[0][rep_el]

        if rep_el.size: 
            atoms.append(int(normal_A[row,7]))
            #print atom_A[row],atom_B[values[0][rep_el]]
            value = np.where(np.add(np.matrix(center_A[row,5],float),np.matrix(center_B[column,5],float)) > np.apply_along_axis(np.linalg.norm,1,np.subtract(np.array(center_A[row,2:5],np.float),np.array(center_B[column,2:5],np.float))))[:][1]

            if value.size:  
                Triangle.append(','.join(normal_A[row]))
                if (condition == 1):
                    temp = atom_A[row]  
                    ligand_atoms = np.concatenate((ligand_atoms,temp.ravel())) 
                #print pdb[atom_A[row]-1],pdb[atom_B[column[value]]+k-1]
                #print atom_A[row],atom_B[column[value]]
                    
                if (temp == atom_A[row]):
                    uniq = np.concatenate((uniq,atom_B[column[value]].ravel()))
                    condition = 0
                    #print temp,uniq
                else :
                    #protein_atoms = np.concatenate((protein_atoms,np.unique(uniq)))
                    condition = 1 
                    uniq = []
                    uniq = np.array(uniq,dtype = int)
                    row = row-1
    #protein_atoms = np.concatenate((protein_atoms,np.unique(uniq)))    
    infoA = []
    infoB = []

    for atom in range(0,len(ligand_atoms)):
        for cenA in range(0,len(center_A)):
            if(int (ligand_atoms[atom]) == int(center_A[cenA,0])):
                infoA.append(center_A[cenA])
                Area = Area + float(center_A[cenA,1])
                break

    #   for cenB in range(0,len(center_B)):     
    #       if(int (protein_atoms[atom]) == int(center_B[cenB,0])):
    #           infoB.append(center_B[cenB])
    #           break
    Proteins = pdb[ligand_atoms-1]
    #print pdb[protein_atoms+k-1]
    infoA = np.array(infoA)
    #distances = np.add(np.matrix(infoA[:,5],float),np.matrix(infoB[:,5],float)) - np.apply_along_axis(np.linalg.norm,1,np.subtract(np.array(infoA[:,2:5],np.float),np.array(infoB[:,2:5],np.float)))
    #depth = np.amax(distances)
    #penetrating_volume(infoA,infoB,distances)
    #print "Penetrating Depth = %f Angstroms" %depth
    return Triangle,Proteins,Area 

def intermediate(normal_A,normal_B,center_A,center_B,j,pdb,Chain):
    if (normal_A.size & normal_B.size):
        Triangle,Proteins,Area = corresponding_molecule(normal_A,normal_B,center_A,center_B,j,pdb,Chain)
        return Triangle,Proteins,Area
    else :return [],0,0


def gepol(l_atom,l_residue,gepol_input,output_path):
    atom_data = []
    normal_data = []
    program_path = os.path.join(os.path.dirname(__file__), 'bin')
    pdbgepol_path = os.path.join(program_path, 'pdbgepol')
    gepol_path = os.path.join(program_path, 'gepol')
    data_file = os.path.join(output_path, 'data.text')
    gepol_file = os.path.join(output_path, 'gepol.text')
    with open(data_file,"w") as myfile:
        start = 1.0 
        gepol_input  = os.path.join(output_path, 'input.text')
        gepol_output = os.path.join(output_path, 'output.text')
        myfile.write("%d\n%d\n%d\n%s\n%s\n%d"%(l_atom,start,l_residue,gepol_input,gepol_output,start)) 

    os.system("{} <{}".format(pdbgepol_path, data_file))
    os.system("{} <{} > {}".format(gepol_path, gepol_output, gepol_file))
    os.remove(gepol_file)
    os.remove(gepol_output)
    os.remove(data_file)
    os.remove("fort.15")
    os.remove(gepol_input)
    with open("fort.7","r") as myfile:
        for _ in xrange(8):
            next(myfile)
        for line in myfile:
            atom_data.append(line.strip())
    os.remove("fort.7")
    with open("fort.8","r") as myfile:
        for _ in xrange(6):
            next(myfile)
        for line in myfile:
            normal_data.append(line.strip())
    os.remove("fort.8")
    return (np.array(atom_data),np.array(normal_data)) 


def formating(info,output_path):
    data = [] 
    length = len(info)
    atoms = 0
    resi_id = info[0][22:26]
    res_number = 1
    gepol_input  = os.path.join(output_path, 'input.text')
    with open(gepol_input,"w") as myfile:
        myfile.write("REMARK PCI\n") 
        for id in range(0,length):
            if (info[id][0:4] == "ATOM"):
                if(resi_id == info[id][22:26]):
                    number = int(res_number)  
                else :
                    res_number = res_number + 1
                    number = int(res_number)
                    resi_id = info[id][22:26] 

                myfile.write("%-9s %-12s %-6d %7s %7s %7s\n"%(info[id][0:6],info[id][12:17],number,info[id][30:38],info[id][38:46],info[id][46:54]))
                line = ("%6s%5s %4s%c%3s %c%4s%c   %8s%8s%8s"%(info[id][0:6],id+1,info[id][12:16],info[id][16],info[id][17:20],info[id][21],number,info[id][26],info[id][30:38],info[id][38:46],info[id][46:54]))
                data.append(line)
                atoms = atoms + 1
            else : continue 
    atom_data,normal_data = gepol(atoms,number,gepol_input,output_path)
    return (np.array(atom_data),np.array(normal_data),np.array(data)) 


def openfile(path): 
    info = []
    with open(path, "r") as fp:
        data = fp.readlines()
        for line in data :
            info.append(line.rstrip())
    return np.array(info)

def statments(A,B,Chain,output_path):
    info_C = []
    info_M = []
    normal = []
    pdb =[]
    atom_A,normal_A,pdb_A  = formating(openfile(A), output_path)
    atom_B,normal_B,pdb_B = formating(openfile(B), output_path)
    molAB = np.concatenate((atom_A,atom_B))
    normal_vector = np.concatenate((normal_A,normal_B))
    pdbAB = np.concatenate((pdb_A,pdb_B))
    atom_AB,normal_AB,pdb_AB = formating(pdbAB, output_path)
    for index in range(0,len(atom_AB)):
        info_C.append(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",atom_AB[index]))
    complexAB = np.array(info_C)[np.where(np.array(info_C)[:,8] == "      0")[:][0]][:]
    for index in range(0,len(molAB)):
        info_M.append(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",molAB[index]))
    moleculeAB = np.array(info_M)[np.where(np.array(info_M)[:,8] == "      0")[:][0]][:]
    for index in range(0,len(normal_vector)):
        normal.append(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",normal_vector[index]))
    normalAB = normal
    for line in pdb_AB:
        pdb.append(line.strip())
    combine_pdb = pdb
    data,i,j = interface(complexAB,moleculeAB)
    normal_A,normal_B,center_A,center_B,j = normal_calculation(data,i,j,normalAB)
    Triangle,Proteins,Area = intermediate(normal_A,normal_B,center_A,center_B,j,pdbAB,Chain)
    #intermediate(normal_B,normal_A,center_B,center_A,0)
    #print np.array(Triangle)
    return Triangle,np.array(Proteins),np.array(Area)

def parse_chain_from_pdb_file_of_one_chain(pdb_file):
    """
    Parse a PDB file of one chain and obtain the atoms for each aminoacid.
    """
    struct = PDB(pdb_file)
    chains = struct.chain_identifiers
    if len(chains) == 1:
        for chain_id in chains:
            chain = struct.get_chain_by_id(chain_id)
    else:
        print('The PDB file {} does not contain one chain!\n'.format(pdb_file))
        sys.exit(10)
    return chain

def renumber_aminoacid_chain_summing_n(chain, n):
    """
    Restart the numeration of the aminoacids summing n to the current number and taking away 1 because it starts from 1.
    """
    for aminoacid in chain.aminoacids:
        aminoacid.number = aminoacid.number + n - 1
    return chain



def main(argv):
    if len(sys.argv)!= 3 :
        print "Usage : python collision_detection.py collision_list_file output_path"
        sys.exit(1)
    else :
        collision_list_file, output_path = argv
        xml_output_file = os.path.join(output_path, 'collision.xml')

        # Read the collision list and execute the calculation for every PPI in the list
        with open(collision_list_file, 'r') as list_fd, open(xml_output_file, 'w') as xml_output_fd:

            xml_output_fd.write('<?xml version="1.0" encoding="utf-8"?>\n<xml>\n')

            for line in list_fd:
                structure_name, A, B = line.strip().split('\t')

                # Replace special characters by other characters (this is also done in the CheckListPPI4.pl script from Baldo)
                structure_name = structure_name.replace('::', '-')
                structure_name = structure_name.replace('#', '_')

                xml_output_fd.write(' <protein>\n')
                xml_output_fd.write('  <id>{}</id>\n'.format(structure_name))
                xml_output_fd.write('  <chain_a>{}</chain_a>\n'.format(os.path.basename(A)))
                xml_output_fd.write('  <chain_b>{}</chain_b>\n'.format(os.path.basename(B)))

                # Parse the chains in the PPI
                chain_A = parse_chain_from_pdb_file_of_one_chain(A)
                chain_B = parse_chain_from_pdb_file_of_one_chain(B)

                Proteins_A = ['ATOM    391  CD  LYS A  51     -97.081  28.268  -7.755','ATOM    392  CE  LYS A  51     -98.178  29.220  -7.318','ATOM    393  NZ  LYS A  51     -98.439  30.125  -8.463','ATOM    678  CG  GLU A  86    -105.261  21.872  -6.521','ATOM    679  CD  GLU A  86    -103.829  21.939  -6.066', 'ATOM   1867  OXT ALA A 232     -93.087  14.468   2.338']
                Proteins_B = ['ATOM    601  O   PRO B  73     -87.598  35.618 -11.563','ATOM    634  SG  CYS B  78     -91.699  33.134  -5.377','ATOM   1319  CG  ASP B 170    -104.317  23.781  -3.722','ATOM   1320  OD1 ASP B 170    -103.057  23.825  -3.728','ATOM   1321  OD2 ASP B 170    -104.987  22.715  -3.714','ATOM   1383  O   LEU B 178     -98.847  31.972  -6.554','ATOM   1468  CG  LEU B 189     -98.788  16.662  -3.629']
                Area_A = '275.181'
                Area_B = '314.459'
                Triangle_A = [12,321]
                Triangle_B = [12,321]

                # Run the calculations
                Triangle_A , Proteins_A , Area_A = statments(A,B,1,output_path)
                Triangle_B , Proteins_B , Area_B = statments(B,A,2,output_path)

                # If the calculations of triangles are not empty, there is collision!
                if Triangle_A != [] and Triangle_B != []:

                    # Output the residues of the proteins
                    residues_A_file = os.path.join(output_path, '{}.residues_A.out'.format(os.path.basename(A)))
                    residues_B_file = os.path.join(output_path, '{}.residues_B.out'.format(os.path.basename(B)))
                    with open(residues_A_file, 'w') as residues_A_fd, open(residues_B_file, 'w') as residues_B_fd:
                        for values in Proteins_A:
                            residues_A_fd.write('{}\n'.format(values))
                        for values in Proteins_B:
                            residues_B_fd.write('{}\n'.format(values))
                    # residues_A_file = "/var/www/html/SPServer_results/ppi_5a153da30c7a8/structA.r_merged.pdb.residues_A.out"
                    # residues_B_file = "/var/www/html/SPServer_results/ppi_5a153da30c7a8/structA.l_merged.pdb.residues_B.out"

                    # Output results in an XML file
                    xml_output_fd.write('  <collision>yes</collision>\n')
                    xml_output_fd.write('  <area_protein_A>{:.3f}</area_protein_A>\n'.format(float(Area_A)))
                    xml_output_fd.write('  <area_protein_B>{:.3f}</area_protein_B>\n'.format(float(Area_B)))

                    for chain_id, chain, residues_file in [[ 'A', chain_A, residues_A_file ], [ 'B', chain_B, residues_B_file ]]:

                        # Parse the aminoacids that collide and output them
                        xml_output_fd.write('  <colliding_residues_{}>\n'.format(chain_id))
                        colliding_residues = parse_chain_from_pdb_file_of_one_chain(residues_file)
                        n = int(chain.first_aminoacid.number)
                        colliding_residues = renumber_aminoacid_chain_summing_n(colliding_residues, n)

                        for residue in colliding_residues.aminoacids:
                            name = residue.type
                            number = residue.number
                            residue_id = '{}-{}'.format(name, str(number))
                            num_atoms_colliding = len(residue)
                            num_atoms_total = len(chain.get_residue_by_identifier(residue.identifier))
                            percentage_atoms_colliding = float(num_atoms_colliding) / float(num_atoms_total) * 100
                            xml_output_fd.write('   <residue>\n')
                            xml_output_fd.write('    <id>{}</id>\n'.format(residue_id))
                            xml_output_fd.write('    <num_atoms_colliding>{}</num_atoms_colliding>\n'.format(num_atoms_colliding))
                            xml_output_fd.write('    <num_atoms_total>{}</num_atoms_total>\n'.format(num_atoms_total))
                            xml_output_fd.write('    <percentage_atoms_colliding>{}</percentage_atoms_colliding>\n'.format(int(round(percentage_atoms_colliding))))
                            xml_output_fd.write('   </residue>\n')
                        xml_output_fd.write('  </colliding_residues_{}>\n'.format(chain_id))

                else:

                    # Output results in an XML file
                    xml_output_fd.write('  <collision>no</collision>\n')
                    xml_output_fd.write('  <area_primary></area_primary>\n')
                    xml_output_fd.write('  <area_secondary></area_secondary>\n')
                    xml_output_fd.write('  <colliding_residues_A></colliding_residues_A>\n')
                    xml_output_fd.write('  <colliding_residues_B></colliding_residues_B>\n')

                xml_output_fd.write(' </protein>\n')

            xml_output_fd.write('</xml>\n')

        print('[SUCCESS] Execution completed!')


if __name__ == "__main__":
   main(sys.argv[1:])

