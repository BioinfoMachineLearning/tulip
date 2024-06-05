"""
processes small molecues because we cannot compute similarity between them. 

"""

import os
from rdkit import Chem
from collections import defaultdict
from openbabel import pybel


p = 'T1158v4' # protein name
small_lig = "MG" #small ligand name
all_lig_dir = f"/protein_ligand/outputs/{p}/" # path where the outputs are saved
protein_na = f"/protein_ligand/may/{p}/query.pdb" # protein 

save_path_mol = f"/protein_ligand/outputs/{p}/{small_lig}" # path where the final ligand outputs are saved
if not os.path.exists(save_path_mol):
    os.makedirs(save_path_mol)
lig_dict = dict()

saved_coords = set()

count = 0

from rdkit.Chem import AllChem
import numpy as np


from Bio.PDB import PDBParser

def extract_coordinates_from_pdb_file_bio(pdb_file, startswith):
    coordinates = list()
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords = atom.get_coord()
                    coordinates.append(coords)
    return coordinates

def extract_coordinates_from_pdb_file(pdb_file, startswith):
    try:
        coordinates = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(startswith):
                    tokens = line.split()
                    x = float(tokens[6])
                    y = float(tokens[7])
                    z = float(tokens[8])
                    coordinates.append([x, y, z])
        return np.array(coordinates)
    except FileNotFoundError:
        return None

def calculate_distance(atom1, atom2):
    return np.linalg.norm(atom1 - atom2)

def is_ligand_on_surface(protein_file, ligand_atoms, threshold):
    try:
        protein_atoms = extract_coordinates_from_pdb_file_bio(protein_file, startswith="ATOM")
        ligand_atoms = [ligand_atoms]
        if protein_atoms is not None and ligand_atoms is not None:
            for ligand_atom in ligand_atoms:
                for protein_atom in protein_atoms:
                    distance = calculate_distance(ligand_atom, protein_atom)
                    if distance <= threshold:
                        return True
    except FileNotFoundError:
        return False

    return False



def read_hetatm_lines(file_path):
    hetatm_lines = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("HETATM"):
                hetatm_lines.append(line.strip())
    return hetatm_lines

def convert_pdb_sdf(PDB_ligand, save_path_mol):
    global count

    lines = read_hetatm_lines(PDB_ligand)
    for line in lines:
        if small_lig in line.split()[2]: 
            l = line.split()[4]
            if len(line.split()[4]) > 1:
                coords_index = 5
            else:
                coords_index = 6
            coords = tuple(map(float, line.split()[coords_index:coords_index+3]))  # Extract coordinates
            if coords not in saved_coords:
                if is_ligand_on_surface(protein_file=protein_na, ligand_atoms=coords, threshold=6):
                    saved_coords.add(coords)
                    cl_data = line.strip()
                    mol = pybel.readstring("pdb", cl_data)
                    
                    # Write the molecule to an SDF file
                    save_path = f"{save_path_mol}/rank_{count}.sdf"
                    mol.write("sdf", save_path, overwrite=True)
                    
                    count += 1
                    print("Done ", count)


    


if __name__ == "__main__":
    all_ligs = pdb_files = [file for file in os.listdir(all_lig_dir) if file.endswith(".pdb")]
    all_ligs.sort()
    for lig in all_ligs:
        li = lig.split("_")
        li = [i for i in li if i!= '']
        lig_name = ''.join(char for char in li[0] if not char.isdigit())
        
        if lig_name == small_lig:
            template_num = li[1].split(".")[0]
            lig_dict[template_num] = lig

    lig_dict = dict(sorted(lig_dict.items()))
    for k, v in lig_dict.items():
        lig_path = f"{all_lig_dir}/{v}"
        convert_pdb_sdf(PDB_ligand=lig_path,save_path_mol=save_path_mol)