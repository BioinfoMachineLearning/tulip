"""
author: nabin 
timestamp: Wed Apr 20 2022 9:56 PM

this runs Equibind to predict the ligand orientations

"""
import shutil
import os
import os.path as osp
import subprocess


def make_equibind_files(input_path, predicted_strus, ligand_file):
    equibind_data = osp.join(input_path, 'equibind_data')
    os.makedirs(equibind_data, exist_ok=True)
    for pred in predicted_strus:
        lig_dir = ligand_file.split(".")[0]
        os.makedirs(osp.join(equibind_data, lig_dir + str(predicted_strus.index(pred))), exist_ok=True)
        shutil.copy(osp.join(input_path, 'predicted', pred),
                    osp.join(equibind_data, lig_dir + str(predicted_strus.index(pred))))
        shutil.copy(osp.join(input_path, ligand_file), osp.join(equibind_data, lig_dir + str(predicted_strus.index(pred))))
    standardize_file(input_path)


def standardize_file(input_path):
    f = [e for e in
         os.listdir(osp.join(input_path, 'equibind_data'))]
    print(f"Equibind files:{f}")
    for i in range(len(f)):
        files_equi = [g for g in os.listdir(osp.join(input_path, 'equibind_data', f[i]))]
        for fs in files_equi:
            if fs.endswith(".sdf"):
                mol_file = fs
                li = fs.split(".")[0]
                li = li.split("_")
                li = li[0] + "_ligand" + ".sdf"
                os.rename(osp.join(input_path, 'equibind_data', f[i], mol_file),
                          osp.join(input_path, 'equibind_data', f[i], li))

            if fs.endswith(".pdb"):
                pdb_file = fs
                pr = fs.split(".")[0]
                pr = pr.split("_")[0]
                pr = pr + "_protein" + ".pdb"
                os.rename(osp.join(input_path, 'equibind_data', f[i], pdb_file),
                          osp.join(input_path, 'equibind_data', f[i], pr))


def run_inference(inference_file, inference_yml):
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Running Equibind")
    proc = subprocess.run(f'python3 {inference_file} --config={inference_yml}', shell=True)
    print('returncode', proc.returncode)
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Equibind Done")
