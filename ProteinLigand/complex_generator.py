"""
author: nabin
timestamp: Wed Apr 20 2022 4:43 PM
"""
import argparse
import os
import os.path as osp

import yaml

from pdbtemplate import check_get_ligand, superimpose, compare_ligand
#from equibindpredictor import predict_ligand_orientation, mol_from_smiles


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=argparse.FileType(mode='r'),
                        default='configs/run.yml')
    parser.add_argument('--input_dir', type=str,
                        default='data/inputs/',
                        help='path where to put the generated complexes')
    parser.add_argument('--output_dir', type=str,
                        default='data/outputs/',
                        help='path where to put the generated complexes')

    return parser.parse_args()


def main():
    global pdb_templates, casp_ligand, pred_stru, lig
    args = parse_arguments()
    if args.config:
        config_dict = yaml.load(args.config, Loader=yaml.FullLoader)
        arg_dict = args.__dict__
        for key, value in config_dict.items():
            if isinstance(value, list):
                for v in value:
                    arg_dict[key].append(v)
            else:
                arg_dict[key] = value
        args.config = args.config.name
    else:
        config_dict = dict()
    print(config_dict['chimera_path'])

    targets = [target for target in os.listdir(config_dict['input_dir'])]
    targets.sort()

    for files in targets:
        if os.path.isdir(osp.join(config_dict['input_dir'], files)):
            predicted_structures = [pred for pred in os.listdir(osp.join(config_dict['input_dir'], files, 'predicted'))]
            casp_ligand = [ligand for ligand in os.listdir(osp.join(config_dict['input_dir'], files)) if
                           ligand.endswith('.smiles')]
            # pdb_templates = [template for template in
                             # os.listdir(osp.join(config_dict['input_dir'], files, 'pdb_templates'))]
            # pdb_templates.sort()
            file = open(osp.join(config_dict['input_dir'], files, 'evalue.m8'))
            pdb_templates = file.readlines()
            casp_ligand.sort()
            predicted_structures.sort()
            casp_ligand = osp.join(config_dict['input_dir'], files, casp_ligand[0])

            print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Found Ligand: {casp_ligand}")
            print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Found PDBTemplate: {file}")
            print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Found Predicted Structures:{predicted_structures}")

            checked_templates = list()

            for tem in range(len(pdb_templates)):
                template = pdb_templates[tem]
                template = template.split("	")[2]
                template = template.split(".")[0].lower()
                template = template[:4]
                template = f"pdb{template}.ent"
                # template_idx = pdb_templates.index(template)
                template_idx = tem
                template = osp.join('templates', template)
                print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Using the template: {template}")
                if template not in checked_templates:
                    checked_templates.append(template)
                    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Checked Templates: ",
                          checked_templates)
                    if check_get_ligand.check_ligands(
                            pdb_template=osp.join(config_dict['input_dir'], files, template)) is True:
                        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Ligand found in template: ",
                              osp.join(config_dict['input_dir'], template))
                        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Running PDBTemplate")

                        for pred_stru in predicted_structures:
                            pred_stru_idx = predicted_structures.index(pred_stru)
                            pred_stru = osp.join('predicted', pred_stru)
                            predicted_stru = osp.join(config_dict['input_dir'],
                                                      files, pred_stru)

                            id_template_pred = f"{template_idx}_{pred_stru_idx}"

                            superimpose.superimpose_structures(chimera_path=config_dict['chimera_path'],
                                                               pdb_template=osp.join(config_dict['input_dir'], files,
                                                                                     template),
                                                               predicted_structure=predicted_stru,
                                                               save_path=osp.join(config_dict['output_dir'], files),
                                                               idx=id_template_pred)
                            check_get_ligand.get_save_ligands(
                                pdb_template=osp.join(config_dict['output_dir'], files,
                                                      f'matchmaker_saved_model_{id_template_pred}.pdb'),
                                save_path=osp.join(config_dict['output_dir'], files,
                                                   f'ligand_template_{id_template_pred}.pdb'))
                            compare_ligand.compute_similarity(PDB_template_ligand=osp.join(config_dict['output_dir'], files,
                                                                                           f'ligand_template_{id_template_pred}.pdb'),
                                                              casp_ligand=casp_ligand, idx=id_template_pred)
                    else:
                        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Ligand not found in template",
                              osp.join(config_dict['input_dir'], template))
                print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Template already done: {template}")

                # running equibind here
                """
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Running Equibind")
                mol_from_smiles.ligands_mol_from_smiles(casp_ligand)
                equibind_data_out = osp.join(config_dict['output_dir'], 'equibind_data_out')
                os.makedirs(equibind_data_out, exist_ok=True)
                lignds = [l for l in os.listdir(osp.join(config_dict['input_dir'], files)) if l.endswith(".mol")]
                print(osp.join(config_dict['input_dir'], files))
                print(f"This is the ligands for equibind {lignds}")
                for lig in lignds:
                    predict_ligand_orientation.make_equibind_files(
                        input_path=osp.join(config_dict['input_dir'], files),
                        predicted_strus=predicted_structures, ligand_file=lig)
                    predict_ligand_orientation.run_inference(inference_file=config_dict['equibind_inference'],
                                                             inference_yml=config_dict['equibind_inference_yml'])
                """

if __name__ == "__main__":
    main()
    print("Process Complete")
