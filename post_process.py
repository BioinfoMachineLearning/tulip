
"""
author: nabin 
timestamp: Tue Oct 17 2023 02:22 PM

- post processing for the extracted ligands. 
- Generate final .sdf file
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter
from openbabel import openbabel
import sys

import os
import subprocess
import time
import re
import numpy as np
from Bio.PDB import PDBParser

def calculate_distance(atom1, atom2):
    return np.linalg.norm(atom1 - atom2)

def extract_coordinates_from_sdf_file(lig_file):
    coordinates = list()
    suppl = Chem.SDMolSupplier(lig_file)
    for mol in suppl:
        if mol:
            conf = mol.GetConformer()
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                coordinates.append((pos.x, pos.y, pos.z))
    return coordinates

def extract_coordinates_from_pdb_file(pdb_file, startswith):
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

def is_ligand_on_surface(protein_file, ligand_file, threshold):
    try:
        protein_atoms = extract_coordinates_from_pdb_file(protein_file, startswith="ATOM")
        ligand_atoms = extract_coordinates_from_sdf_file(ligand_file)
        if protein_atoms is not None and ligand_atoms is not None:
            for ligand_atom in ligand_atoms:
                for protein_atom in protein_atoms:
                    distance = calculate_distance(ligand_atom, protein_atom)
                    if distance <= threshold:
                        return True
    except FileNotFoundError:
        return False

    return False

total_rank_ligands_to_return = 49

def extract_query_mol(output_dir):
    pdb_files = [file for file in os.listdir(output_dir) if file.endswith(".pdb")]
    pdb_files.sort()
    count_q = 1
    done_files = list()
    for f in pdb_files:
        fi = f"{output_dir}/{f}"
        with open(fi, 'r') as file:
            contents = file.readlines()
        current_section = list()
        for line in contents:
            current_section.append(line)
            if "TER\n" in line:
                filename = f'{output_dir}/rank_{count_q}.pdb'
                done_files.append(filename)
                with open(filename, 'w') as output_file:
                    output_file.writelines(current_section)
                output_file.close()
                count_q += 1
                break
    return done_files

def pdb_mol(pdb_file, mol_file):
    mol = Chem.MolFromPDBFile(pdb_file)
    if mol is not None:
        Chem.MolToMolFile(mol, mol_file)

def mol_mol2(mol_file, mol2_file):
    obabel_path = "/bml/nabin/miniconda_li/envs/tulip/bin/obabel"
    cmd = [obabel_path, mol_file, "-O", mol2_file]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    stdout = result.stdout
    stderr = result.stderr
    return_code = result.returncode
    if return_code == 0:
        print("Openbabel executed successfully.")
        print(stdout)
    else:
        print(f"Openbabel failed with exit code {return_code}.")
        print("Standard Error:")
        print(stderr)

def save_mol(query_smiles, query_mol_save_path):
    mol = Chem.MolFromSmiles(query_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)  # Adjust the seed
    Chem.MolToMolFile(mol, query_mol_save_path)


def generate_sdf(output_dir, mol_files):
    count = 1
    for i in mol_files:
        mol_file_path = i
        sdf_file_path = f"{output_dir}/rank_{count}.sdf"
        if os.path.exists(sdf_file_path):
            os.remove(sdf_file_path)
        try:
            mol_supplier = Chem.SDMolSupplier(mol_file_path)
            sdf_writer = SDWriter(sdf_file_path)
            for idx, mol in enumerate(mol_supplier):
                if mol is not None:
                    mol.SetProp("_Name", f"Molecule_{idx + 1}")
                    sdf_writer.write(mol)
                    count += 1
            sdf_writer.close()
        except OSError:
            pass


def rank_them(output_dir, query_protein):
    mol_sim = dict()
    all_files = os.listdir(output_dir)

    for f in all_files:
        try:
            mol_similarity = f.split("_")[2]
            mol_sim[f] = mol_similarity
        except IndexError:
            pass
    if len(mol_sim) != 0:
        sorted_mol_sim_dict = dict(sorted(mol_sim.items(), key=lambda item: item[1], reverse=True))  
        top_n = {k: sorted_mol_sim_dict[k] for k in list(sorted_mol_sim_dict)[:total_rank_ligands_to_return]}
        return top_n
    else:
        return None


def remove_intermediatory_files(output_dir, retain_extension):
    all_files = os.listdir(output_dir)
    filtered_files = [file for file in all_files if file.endswith(retain_extension)]

    # Remove files with other extensions
    for file in all_files:
        if file not in filtered_files:
            file_path = os.path.join(output_dir, file)
            os.remove(file_path)  # Remove the file



def run_LSalign(query_smiles, output_dir, top_n, ls_align_path):
    """
    query_ligand = mol file
    template_ligand = mol file
    """
    

    query_mol2_save_path = f"{output_dir}/query_ligand.mol2"
    query_mol_save_path = f"{output_dir}/query_ligand.mol"
    save_mol(query_smiles, query_mol_save_path)
    mol_mol2(mol_file=query_mol_save_path, mol2_file=query_mol2_save_path)


    count = 1
    for key in top_n:
        mol_file_path = f"{output_dir}/{key}"
        template_mol2_save_path = f"{output_dir}/{count}_template_ligand.mol2"
        mol_mol2(mol_file=mol_file_path, mol2_file=template_mol2_save_path)

        time.sleep(5)
        
        output_filename = f"{output_dir}/{count}_aligned.pdb"
        cmd = [ls_align_path, query_mol2_save_path, template_mol2_save_path, "-o", output_filename, "-rf", "1", "-acc", "1", "-md", "1" ]
        print(cmd)
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            stdout = result.stdout
            stderr = result.stderr
            return_code = result.returncode
            count += 1
            if return_code == 0:
                print("LS align executed successfully.")
                print(stdout)
            else:
                print(f"LS align failed with exit code {return_code}.")
                print("Standard Error:")
                
                print(stderr)

        except subprocess.CalledProcessError as e:
            # print(f"Command failed with exit code {e.returncode}.")
            print(e.stderr)
            count += 1
        except Exception as e:
            print(f"An error occurred: {str(e)}")


if __name__ == "__main__":
    protein_id = sys.argv[1]

    # protein_id = "1G9V_RQ3"

    ligand_id = "0" # ligand id 

    
    output_dir1 = f"/output/{protein_id}/{ligand_id}" # ligand id dir path
    query_ligand_file = f"/ligands/{protein_id}_unique.smiles" # ligand path

    query_protein = f"/templates/{protein_id}/query.pdb" # query ligand

    print("Working:", protein_id, ligand_id)
    print(query_ligand_file)

    query_ligand = None
    with open(query_ligand_file, 'r') as file:
        for line in file:
            compound_data = re.split(r'\s+', line.strip())
            compound_data = [data for data in compound_data if data]
            smiles_string, ID = compound_data[1], compound_data[0]
            if ID == ligand_id:
                query_ligand = smiles_string
            


    
    ls_align_path = "/bml/tools/LSalign/src/LSalign"
    top_n = rank_them(output_dir1, query_protein)
    remove_intermediatory_files(output_dir1, retain_extension=".mol")

    if top_n != None:
        run_LSalign(query_ligand, output_dir1, top_n, ls_align_path)
        done_files = extract_query_mol(output_dir1)
        done_files.sort()
        mol_files = list()
        for i in range(len(done_files)):
            mol_file = f"{output_dir1}/rank_{i + 1}.mol"
            pdb_mol(done_files[i], mol_file)
            mol_files.append(mol_file)
        generate_sdf(output_dir1, mol_files)
        remove_intermediatory_files(output_dir1, retain_extension=".sdf")
        print("Done", i)
    else: 
        print("=>>>>>>>>>>>>>>>>>>>>>>No ligands found! =>>>>>>>>>>>>>>>>>>>>>>> ",i)
    print("DONE ALL")        
