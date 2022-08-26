"""
author: nabin 
timestamp: Tue May 03 2022 9:22 PM

- this compares similarity between ligands in SMILES format
"""
import os

from pyrosetta import *

init()

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from pdbtemplate import save_pdb_format as save_pdb

similarity_threshold = 0.0


def read_smiles(casp_ligand):
    smiles_found = list()
    lig_names = set()
    ligand_file = open(casp_ligand)
    ligand_smiles = ligand_file.readlines()
    for line in ligand_smiles:
        line_split = line.split(" ")
        print(line_split)
        print(len(line_split))
        # print(line_split[0])
        # print(line_split[1])
        # print(line_split[2])
        lig_smi = line_split[2]
        # lig_smi = lig_smi.split(" ")
        print(lig_smi)
        lig_name = line_split[1]
        smi = lig_smi
        print(f"CASP Ligand Name: {lig_name}")
        # if lig_name == "LIG":
        # smiles_found.append(smi.strip("\n"))
        smiles_found.append(smi)
        # lig_names.add(lig_name.strip("\n"))
        lig_names.add(lig_name)
    #print(smiles_found)
    print(lig_names)
    return smiles_found, lig_names


def save_ligands_seperate(PDB_template_ligand, ligands_found, idx):
    # os.makedirs(save_path, exist_ok=True)
    for ligs in ligands_found:
        parent = os.path.join(PDB_template_ligand, os.pardir)
        save_path = os.path.join(os.path.abspath(parent),
                                 f"{ligs}_{idx}.pdb")
        count = 0
        print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Extracting Ligands {ligs}")
        pose_pdb = pose_from_pdb(PDB_template_ligand)
        for i in range(1, pose_pdb.total_residue() + 1):
            if pose_pdb.residue(i).is_ligand():
                if pose_pdb.residue(i).name()[-3:] == ligs:
                    print(f"The identified ligand is :{pose_pdb.residue(i).name()}")
                    print(f"The identified chain of the ligand is :{pose_pdb.pdb_info().chain(i)}")
                    for atm in range(1, pose_pdb.residue(i).natoms() + 1):
                        x = round(pose_pdb.residue(i).xyz(atm)[0], 3)
                        y = round(pose_pdb.residue(i).xyz(atm)[1], 3)
                        z = round(pose_pdb.residue(i).xyz(atm)[2], 3)
                        atom_name = pose_pdb.residue(i).atom_name(atm).strip()
                        ligand_name = pose_pdb.residue(i).name()[-3:]
                        chain_name = pose_pdb.pdb_info().chain(i)
                        num_chain = pose_pdb.pdb_info().pose2pdb(i)
                        res_num = int(num_chain.split(" ")[0])
                        save_pdb.save(file_path=save_path, x_cord=x, y_cord=y, z_cord=z, atom_name=atom_name,
                                      ligand_name=ligand_name, chain_id=chain_name, residue_number=res_num, count=count)
                        count += 1
                    save_pdb.save_TER(file_path=save_path)
        save_pdb.save_end(file_path=save_path)


def check_multiple_ligands(PDB_template_ligand, idx):
    ligands_found = set()
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Checking for Multiple Ligands in PDB template Ligand File")
    pose_pdb = pose_from_pdb(PDB_template_ligand)
    # print(f"Total number of residue found for the pdb protein is : {pose_pdb.total_residue()}")
    for i in range(1, pose_pdb.total_residue() + 1):
        if pose_pdb.residue(i).is_ligand():
            lig_name = pose_pdb.residue(i).name().strip()
            ligands_found.add(lig_name[-3:])
    print(
        f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Total Ligands Identified in Multi Ligand Checks are: {ligands_found}")
    save_ligands_seperate(PDB_template_ligand, ligands_found, idx)
    return ligands_found


def compute_similarity(PDB_template_ligand, casp_ligand, idx):
    """

    :param PDB_template_ligand: pdb file saved relative to predicted structure
    :param casp_ligand: smiles file from casp
    :return:
    """
    smiles_list, lig_names = read_smiles(casp_ligand)
    if os.path.exists(PDB_template_ligand):
        parent = os.path.join(PDB_template_ligand, os.pardir)
        lig_names = check_multiple_ligands(PDB_template_ligand, idx)

        for ligs in lig_names:
            template_ligand = os.path.join(os.path.abspath(parent),
                                           f"{ligs}_{idx}.pdb")
            PDB_lig_mol = Chem.MolFromPDBFile(molFileName=template_ligand, removeHs=True)
            # PDB_mol_smiles = Chem.MolToSmiles(PDB_lig_mol)
            try:
                ref_ECFP4_fps = AllChem.GetMorganFingerprintAsBitVect(PDB_lig_mol, 2)
                for smi in smiles_list:
                    print(smi)
                    check_mol = Chem.MolFromSmiles(smi)
                    check_ECFP4_fps = AllChem.GetMorganFingerprintAsBitVect(check_mol, 2)
                    similarity = DataStructs.FingerprintSimilarity(ref_ECFP4_fps,
                                                                   check_ECFP4_fps) # default FingerprintSimilarity is Tanimoto similarity
                    # similarity1 = DataStructs.TverskySimilarity(ref_ECFP4_fps, check_ECFP4_fps)
                    similarity2 = DataStructs.DiceSimilarity(ref_ECFP4_fps, check_ECFP4_fps)
                    # similarity3 = DataStructs.RusselSimilarity(ref_ECFP4_fps, check_ECFP4_fps)
                    # similarity4 = DataStructs.CosineSimilarity(ref_ECFP4_fps, check_ECFP4_fps)
                    print(
                        f">>>>>>>>>>> Similarity between PDB Template Ligand {template_ligand} and CASP Ligand {smiles_list.index(smi)} is: {round(similarity, 3)} <<<<<<<<<<<<<")
                    if similarity >= similarity_threshold:
                        parent = os.path.join(PDB_template_ligand, os.pardir)
                        save_path = os.path.join(os.path.abspath(parent),
                                                 f"{ligs}_mol_sim_Finger_{round(similarity, 3)}_Dice_{round(similarity2, 3)}_{idx}_smiles_{smiles_list.index(smi)}.mol")
                        print(Chem.MolToMolBlock(PDB_lig_mol),
                              file=open(save_path, 'w+'))
                    else:
                        print(
                            f">>>>>>>>>>> Similarity {round(similarity, 3)} for template_predicted {idx} and CASP ligand {smiles_list.index(smi)} is less than threshold: {similarity_threshold} <<<<<<<<<<<<<")

            except (RuntimeError, TypeError, NameError, BaseException) as Error:
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("Fingerprint ran into error for", template_ligand)
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

                continue
