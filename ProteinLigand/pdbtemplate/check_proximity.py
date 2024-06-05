"""
author: nabin 
timestamp: Tue Oct 16 2023 11:22 AM

- checks if ligand is in protein surface. returns bool : True or False
"""


from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


from Bio.PDB import PDBParser


def extract_coordinates_from_pdb_file_bio(pdb_file, startswith):
    coordinates = list()
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)

    # Iterate over all atoms in the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # Get the atom's coordinates
                    coords = atom.get_coord()
                    # Append the coordinates to the list
                    coordinates.append(coords)
    return coordinates

def extract_coordinates_from_pdb_file(pdb_file, startswith):
    try:
        coordinates = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(startswith):
                    tokens = line.split()
                    try:
                        x = float(tokens[5])
                        y = float(tokens[6])
                        z = float(tokens[7])
                        coordinates.append([x, y, z])
                    except ValueError:
                        pass
        return np.array(coordinates)
    except FileNotFoundError:
        return None

def get_coordinates_from_pdb(pdb_file):
    try:
        molecule = Chem.MolFromPDBFile(pdb_file)
        if molecule is not None:
            conf = molecule.GetConformer(0)
            return np.array([conf.GetAtomPosition(i) for i in range(molecule.GetNumAtoms())])
        return None
    except OSError:
        return None

def calculate_distance(atom1, atom2):
    return np.linalg.norm(atom1 - atom2)

def is_ligand_on_surface(protein_file, ligand_file, threshold):
    try:
        protein_atoms = extract_coordinates_from_pdb_file_bio(protein_file, startswith="ATOM")
        ligand_atoms = extract_coordinates_from_pdb_file(ligand_file, startswith="HETATM")
        if protein_atoms is not None and ligand_atoms is not None:
            for ligand_atom in ligand_atoms:
                for protein_atom in protein_atoms:
                    distance = calculate_distance(ligand_atom, protein_atom)
                    if distance <= threshold:
                        return True
    except FileNotFoundError:
        return False

    return False

