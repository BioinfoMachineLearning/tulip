"""
author: nabin 
timestamp: Tue May 03 2022 9:22 PM

- this compares similarity between ligands in SMILES format
"""
import os


from pyrosetta import init, pose_from_pdb
init(options="-mute all")

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from pdbtemplate import save_pdb_format as save_pdb
from pdbtemplate import check_proximity, separate_identical_ligands


similarity_threshold = 0.1
protein_ligand_threshold = 6.0


def read_smiles(casp_ligand):
    smiles_found = dict()
    ligand_file = open(casp_ligand)
    ligand_smiles = ligand_file.readlines()
    for line in ligand_smiles:
        line_split = line.split("\t")
        smiles_found[line_split[1].strip("\n")] = line_split[0].strip("\t")
    return smiles_found


def save_ligands_seperate(PDB_template_ligand, ligands_found, idx):
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
                        atom_name = pose_pdb.residue(i).atom_name(atm)
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
    try:
        pose_pdb = pose_from_pdb(PDB_template_ligand)
        print(f"Total number of residue found for the pdb protein is : {pose_pdb.total_residue()}")
        for i in range(1, pose_pdb.total_residue() + 1):
            if pose_pdb.residue(i).is_ligand():
                ligands_found.add(pose_pdb.residue(i).name()[-3:])
        print(
            f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Total Ligands Identified in Multi Ligand Checks are: {ligands_found}")
        save_ligands_seperate(PDB_template_ligand, ligands_found, idx)
    except RuntimeError as e:
        return ligands_found
    return ligands_found


def compute_similarity(PDB_template_ligand, casp_ligand, query_protein, idx):
    """

    :param PDB_template_ligand: pdb file saved relative to predicted structure
    :param casp_ligand: smiles file from casp
    :return:
    """
    lig_names1 = check_multiple_ligands(PDB_template_ligand, idx)
    lig_names = list()
    for l in lig_names1:
        parent = os.path.join(PDB_template_ligand, os.pardir)
        template_ligand = os.path.join(os.path.abspath(parent),f"{l}_{idx}.pdb")
        lig_names2 = separate_identical_ligands.sep_identical_ligs(template_ligand, idx)
        lig_names.extend(lig_names2)

    for ligs in lig_names:
        parent = os.path.join(PDB_template_ligand, os.pardir)
        template_ligand = os.path.join(os.path.abspath(parent),
                                       f"{ligs}_{idx}.pdb")


        is_on_surface = check_proximity.is_ligand_on_surface(protein_file=query_protein, ligand_file=template_ligand, threshold=protein_ligand_threshold)
        if is_on_surface:
            PDB_lig_mol = Chem.MolFromPDBFile(molFileName=template_ligand, removeHs=True)
            if PDB_lig_mol is not None:
                ref_ECFP4_fps = AllChem.GetMorganFingerprintAsBitVect(PDB_lig_mol, 2)
                smiles_d = read_smiles(casp_ligand)
                for smi, ID in smiles_d.items():
                    check_mol = Chem.MolFromSmiles(smi)
                    check_ECFP4_fps = AllChem.GetMorganFingerprintAsBitVect(check_mol, 2)
                    similarity = DataStructs.FingerprintSimilarity(ref_ECFP4_fps, check_ECFP4_fps)
                    print(
                        f"Similarity between PDB Template Ligand {template_ligand} and CASP Ligand {ID} is: {round(similarity,3)}")
                    if similarity >= similarity_threshold:
                        parent = os.path.join(PDB_template_ligand, os.pardir)
                        save_path_dir = os.path.join(os.path.abspath(parent), ID)
                        if not os.path.exists(save_path_dir):
                            os.makedirs(save_path_dir)
                        save_path = os.path.join(save_path_dir, 
                                                f"{ligs}_sim_{round(similarity,3)}_{idx}.mol")
                        print(Chem.MolToMolBlock(PDB_lig_mol),
                            file=open(save_path, 'w+'))
                    else:
                        print(
                            f"Similarity {round(similarity,3)} for template_predicted {idx} and CASP ligand {ID} is less than threshold: {similarity_threshold}")
