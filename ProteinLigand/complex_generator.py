"""
author: nabin
timestamp: Wed Apr 20 2022 4:43 PM
Updated : Jan 2024


"""
import argparse
import os
import os.path as osp
import gzip
import warnings
from rdkit import Chem
import pandas as pd
import re

# Suppress all warnings
warnings.filterwarnings("ignore")


import yaml
COMMENT_MARKER = '#'
from pdbtemplate import check_get_ligand, superimpose, compare_ligand



def process_templates(template_path, template_hits, query_protein, query_ligand, chimera_path, output_dir):

    total_templates = len(template_hits)
    count = 0

    for idx in range(len(template_hits)):
        template = f"{template_path}/{template_hits[idx]}"
        if check_get_ligand.check_ligands(pdb_template=template) is True:
            output_name = os.path.basename(query_protein).split(".")[0]
            output_path = f"{output_dir}/{output_name}"
            superimpose.superimpose_structures(chimera_path=chimera_path, pdb_template=template, predicted_structure=query_protein, 
            save_path=output_path, idx=idx)
            matchmaker_savename = osp.join(output_path, f'matchmaker_saved_model_{idx}.pdb')
            ligand_template_savename = osp.join(output_path, f'ligand_template_{idx}.pdb')
            check_get_ligand.get_save_ligands(pdb_template=matchmaker_savename, save_path=ligand_template_savename)
            compare_ligand.compute_similarity(PDB_template_ligand=ligand_template_savename, casp_ligand=query_ligand, query_protein=query_protein, idx=idx)
            print("Extracted Ligand from template :", template_hits[idx])
            count += 1
            print(f"Completed scanning templates: {count} / {total_templates}")

            os.remove(matchmaker_savename)
            os.remove(ligand_template_savename)
            
        else:
            count += 1
            print(f"Completed scanning templates: {count} / {total_templates}")
            

def unzip_files(protein_dir, templates):
    for gz_file in templates:
        gz_path = f"{protein_dir}/{gz_file}"
        file_name = gz_file.split(".")[0]
        unzipped_file = f"{protein_dir}/{file_name}.pdb" 
        with gzip.open(gz_path, 'rb') as gz_file:
            with open(unzipped_file, 'wb') as output_file:
                output_file.write(gz_file.read())
        os.remove(gz_path)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=argparse.FileType(mode='r'),
                        default='tulip/configs/configs.yml')
    parser.add_argument('--protein_name', type=str)
    parser.add_argument('--output_dir', type=str,
                        default='outputs',
                        help='path where to put the generated complexes')


    return parser.parse_args()


def parse_args(args):
    global templates_hits, ligand_smiles, query_structure, lig
    if args.config is not None:
        config_dict = yaml.safe_load(args.config)
        config_dict = {k: v for k, v in config_dict.items() if not k.startswith(COMMENT_MARKER)}
        args.config = args.config.name
    else:
        config_dict = dict()
    
    if args.protein_name is not None:
        config_dict['protein_name'] = args.protein_name

    print(config_dict)

    
    if args.protein_name is not None:
        config_dict['protein_name'] = args.protein_name

    print(config_dict)
    print(args.protein_name)
    template_dir = config_dict['template_dir']
    protein_templates = config_dict['protein_templates']
    output_dir = config_dict['output_dir']
    chimera_path = config_dict['chimera_path']
    protein_name = config_dict['protein_name']
    ligand_dir = config_dict['ligand_dir']

    return template_dir, protein_templates, output_dir, chimera_path, protein_name, ligand_dir


def main(template_dir, protein_templates, output_dir, chimera_path, protein, ligand_dir):
    query_protein = f"{protein_templates}/{protein}/{protein}.pdb"
    query_ligands_file = f"{ligand_dir}/{protein}.smiles.txt"
    template_names_csv_file = f"{protein_templates}/{protein}/template_names.csv"
    smiles_dict = dict()
    with open(query_ligands_file, 'r') as file:
        next(file)
        for line in file:
            compound_data = re.split(r'\s+', line.strip())
            compound_data = [data for data in compound_data if data]
            smiles_string, ID = compound_data[2], compound_data[0]
            # print(compound_data)

            if smiles_string not in smiles_dict:
                smiles_dict[smiles_string] = ID
    print(smiles_dict)
    print("Unique Ligands: ", len(smiles_dict))
    unique_ligand_file = f'{ligand_dir}/{protein}_unique.smiles'
    with open(unique_ligand_file, 'w') as s:
        for k, v in smiles_dict.items():
            s.write(f'{v}\t{k}')
            s.write('\n')

    template_names_df = pd.read_csv(template_names_csv_file, delimiter='\t')
    df_sorted = template_names_df.sort_values(by='evalue', ascending=True)
    template_hits = {ele[:4].lower() + ".pdb1" for ele in df_sorted['target']}
    template_hits = list(template_hits)
    print("##############################")
    print("Working with: ")
    print("Ligand file :", unique_ligand_file)
    print("Query structure :", query_protein)
    print("Total template hits :", len(template_hits))
    print("##############################")
    
    process_templates(template_dir, template_hits, query_protein, unique_ligand_file, chimera_path, output_dir)

if __name__ == "__main__":
    args = parse_arguments()
    template_dir, protein_templates, output_dir, chimera_path, protein_name, ligand_dir = parse_args(args)
    main(template_dir, protein_templates, output_dir, chimera_path, protein_name, ligand_dir)
    print("Process Complete")
