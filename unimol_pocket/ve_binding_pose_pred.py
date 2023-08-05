"""
main file for binding pose prediction
code adapted from Uni-Mol Docking Colab
https://colab.research.google.com/github/dptech-corp/Uni-Mol/blob/main/unimol/notebooks/unimol_binding_pose_demo.ipynb
"""

import os
import numpy as np
import pandas as pd
import lmdb
from biopandas.pdb import PandasPdb
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.cluster import KMeans
from rdkit.Chem.rdMolAlign import AlignMolConformers
from unimol.utils.docking_utils import docking_data_pre, ensemble_iterations
import pickle
import re
import json
import copy
import concurrent.futures
import subprocess


main_atoms = ["N", "CA", "C", "O", "H"]

def load_from_data(pdb_id, input_dir): # this is pdbid_pid
    try:
        pdbid, _ = pdb_id.split("_")
        pdb_path = os.path.join(input_dir, "veprotein", pdbid + ".first_chain_out.pdb")
        pmol = PandasPdb().read_pdb(pdb_path)
        pocket_residues = json.load(
            open(os.path.join(input_dir, "veprotein.pocket.json"))
        )[pdb_id]
        return pmol, pocket_residues
    except:
        print("Currently not support parsing pdb and pocket info from local files.")


def normalize_atoms(atom):
    return re.sub("\d+", "", atom)


def single_conf_gen(tgt_mol, num_confs=1000, seed=42, removeHs=True):
    mol = copy.deepcopy(tgt_mol)
    mol = Chem.AddHs(mol)
    allconformers = AllChem.EmbedMultipleConfs(
        mol, numConfs=num_confs, randomSeed=seed, clearConfs=True
    )
    sz = len(allconformers)
    for i in range(sz):
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=i)
        except:
            continue
    if removeHs:
        mol = Chem.RemoveHs(mol)
    return mol


def clustering_coords(mol, M=1000, N=100, seed=42, removeHs=True):
    rdkit_coords_list = []
    rdkit_mol = single_conf_gen(mol, num_confs=M, seed=seed, removeHs=removeHs)
    noHsIds = [
        rdkit_mol.GetAtoms()[i].GetIdx()
        for i in range(len(rdkit_mol.GetAtoms()))
        if rdkit_mol.GetAtoms()[i].GetAtomicNum() != 1
    ]
    ### exclude hydrogens for aligning
    AlignMolConformers(rdkit_mol, atomIds=noHsIds)
    sz = len(rdkit_mol.GetConformers())
    for i in range(sz):
        _coords = rdkit_mol.GetConformers()[i].GetPositions().astype(np.float32)
        rdkit_coords_list.append(_coords)

    ### exclude hydrogens for clustering
    rdkit_coords_flatten = np.array(rdkit_coords_list)[:, noHsIds].reshape(sz, -1)
    ids = (
        KMeans(n_clusters=N, random_state=seed)
        .fit_predict(rdkit_coords_flatten)
        .tolist()
    )
    coords_list = [rdkit_coords_list[ids.index(i)] for i in range(N)]
    return coords_list


def parser(pdb_id, smiles, input_dir, seed=42):
    pmol, pocket_residues = load_from_data(pdb_id, input_dir)
    pname = pdb_id
    pro_atom = pmol.df["ATOM"]
    pro_hetatm = pmol.df["HETATM"]

    pro_atom["ID"] = pro_atom["chain_id"].astype(str) + pro_atom[
        "residue_number"
    ].astype(str)
    pro_hetatm["ID"] = pro_hetatm["chain_id"].astype(str) + pro_hetatm[
        "residue_number"
    ].astype(str)

    pocket = pd.concat(
        [
            pro_atom[pro_atom["ID"].isin(pocket_residues)],
            pro_hetatm[pro_hetatm["ID"].isin(pocket_residues)],
        ],
        axis=0,
        ignore_index=True,
    )

    pocket["normalize_atom"] = pocket["atom_name"].map(normalize_atoms)
    pocket = pocket[pocket["normalize_atom"] != ""]
    patoms = pocket["atom_name"].apply(normalize_atoms).values.tolist()
    pcoords = [pocket[["x_coord", "y_coord", "z_coord"]].values]
    side = [0 if a in main_atoms else 1 for a in patoms]
    residues = (
        pocket["chain_id"].astype(str) + pocket["residue_number"].astype(str)
    ).values.tolist()

    M, N = 100, 10
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    latoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    holo_coordinates = [mol.GetConformer().GetPositions().astype(np.float32)]
    holo_mol = mol
    coordinate_list = clustering_coords(mol, M=M, N=N, seed=seed, removeHs=False)
    mol_list = [mol] * N

    return pickle.dumps(
        {
            "atoms": latoms,
            "coordinates": coordinate_list,
            "mol_list": mol_list,
            "pocket_atoms": patoms,
            "pocket_coordinates": pcoords,
            "side": side,
            "residue": residues,
            "holo_coordinates": holo_coordinates,
            "holo_mol": holo_mol,
            "holo_pocket_coordinates": pcoords,
            "smi": smiles,
            "pocket": pname,
        },
        protocol=-1,
    )


def write_lmdb(pdb_id, smiles_list, input_dir, seed=42, result_dir="./results"):
    os.makedirs(result_dir, exist_ok=True)
    outputfilename = os.path.join(result_dir, pdb_id + ".lmdb")
    try:
        os.remove(outputfilename)
    except:
        pass
    env_new = lmdb.open(
        outputfilename,
        subdir=False,
        readonly=False,
        lock=False,
        readahead=False,
        meminit=False,
        max_readers=1,
        map_size=int(10e9),
    )
    for i, smiles in enumerate(smiles_list):
        inner_output = parser(pdb_id, smiles, input_dir, seed=seed)
        txn_write = env_new.begin(write=True)
        txn_write.put(f"{i}".encode("ascii"), inner_output)
    txn_write.commit()
    env_new.close()


def generate_docking_input(
    predict_file, reference_file, tta_times=10, output_dir="./results"
):
    (
        mol_list,
        smi_list,
        pocket_list,
        pocket_coords_list,
        distance_predict_list,
        holo_distance_predict_list,
        holo_coords_list,
        holo_center_coords_list,
    ) = docking_data_pre(reference_file, predict_file)
    iter = ensemble_iterations(
        mol_list,
        smi_list,
        pocket_list,
        pocket_coords_list,
        distance_predict_list,
        holo_distance_predict_list,
        holo_coords_list,
        holo_center_coords_list,
        tta_times=tta_times,
    )
    for i, content in enumerate(iter):
        pocket = content[3]
        output_name = os.path.join(output_dir, "{}.{}.pkl".format(pocket, i))
        try:
            os.remove(output_name)
        except:
            pass
        pd.to_pickle(content, output_name)


def generate_binding_pose(smiles_list, data_path_meta, weight_path, seed=42, batch_size=1, dist_threshold=8.0, recycling=3):
    
    """
    smiles_list: list of smiles
    data_path_meta: input data dir of `generate_binding_pose`
    results_path: output dir of `generate_binding_pose`
    weight_path: path of pretrained model

    Note that all `pdb_id` here are actually pdbid_pocketid
    """

    data_path = os.path.join(data_path_meta, "protein.lmdb/")
    results_path = os.path.join(data_path_meta, "results/")
    pdbid_list = list(json.load(open(os.path.join(data_path_meta, "veprotein.pocket.json"))).keys())

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(write_lmdb, pdb_id, smiles_list, data_path_meta, seed=seed, result_dir=data_path) for pdb_id in pdbid_list]
        for future, item in zip(concurrent.futures.as_completed(futures), pdbid_list):
            try:
                result = future.result()
                print(item)
            except Exception as e:
                print(f"Error while processing item {item}: {str(e)}")

    def run_bash_command(command):
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        return result.stdout, result.stderr, result.returncode

    commands = []
    for pdb_id in json.load(open(os.path.join(data_path_meta, "veprotein.pocket.json"))).keys():
        device = np.random.choice([0, 1, 5, 6, 7])
        commands.append(f"CUDA_VISIBLE_DEVICES={device} python ./unimol/infer.py --user-dir ./unimol {data_path} --valid-subset {pdb_id} \
        --results-path {results_path} \
        --num-workers 8 --ddp-backend=c10d --batch-size {batch_size} \
        --task docking_pose --loss docking_pose --arch docking_pose \
        --path {weight_path} \
        --fp16 --fp16-init-scale 4 --fp16-scale-window 256 \
        --dist-threshold {dist_threshold} --recycling {recycling} \
        --log-interval 50 --log-format simple")

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(run_bash_command, cmd) for cmd in commands]
        for future, command in zip(concurrent.futures.as_completed(futures), commands):
            try:
                stdout, stderr, returncode = future.result()
                if returncode == 0:
                    print(f"Command '{command}' executed successfully.")
                    print("STDOUT:", stdout)
                else:
                    print(f"Command '{command}' failed with error:")
                    print("STDOUT:", stdout)
                    print("STDERR:", stderr)
            except Exception as e:
                print(f"Error while executing command '{command}': {str(e)}")

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for pdb_id in pdbid_list:
            predict_file = os.path.join(results_path, ".._" + pdb_id + ".out.pkl")
            reference_file = os.path.join(data_path, pdb_id + ".lmdb")
            futures.append(executor.submit(generate_docking_input, predict_file, reference_file, tta_times=10, output_dir=results_path))
        for future, item in zip(concurrent.futures.as_completed(futures), pdbid_list):
            try:
                result = future.result()
                print(item)
            except Exception as e:
                print(f"Error while processing item {item}: {str(e)}")

    docking_commands = []
    for pdb_id in pdbid_list:
        input_path = os.path.join(results_path, "{}.0.pkl".format(pdb_id))
        ligand_path = os.path.join(results_path, "docking.{}.sdf".format(pdb_id))
        docking_commands.append(f"python ./unimol/utils/coordinate_model.py --input {input_path} --output-ligand {ligand_path}")

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(run_bash_command, cmd) for cmd in docking_commands]
        for future, command in zip(concurrent.futures.as_completed(futures), docking_commands):
            try:
                stdout, stderr, returncode = future.result()
                if returncode == 0:
                    print(f"Command '{command}' executed successfully.")
                    print("STDOUT:", stdout)
                else:
                    print(f"Command '{command}' failed with error:")
                    print("STDOUT:", stdout)
                    print("STDERR:", stderr)
            except Exception as e:
                print(f"Error while executing command '{command}': {str(e)}")
