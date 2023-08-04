from pathlib import Path
import numpy as np
import json
import os
from utils import run_command, set_directory
from unimol_pocket.op import tleap, create_system, sdf2pdb
from amber import em


def calc_affinity(mmpbsa_in, cp, rp, lp, crd):
    """Calculate the binding affinity."""
    # MMPBSA.py -O -i mmpbsa.in -cp com.top -rp rec.top -lp lig.top -y traj.crd
    run_command(
        ["MMPBSA.py", "-O", "-i", mmpbsa_in, "-cp", cp, "-rp", rp, "-lp", lp, "-y", crd]
    )
    return

pocket_json = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data/VEprotein/veprotein.pocket.json")
result_dir = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data/VEprotein/results")
affinity_dir = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data/VEprotein/affinity")
ligand_topology = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data/VEprotein/ligand_topology")
protein_processed_pdbs = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data/processed")

previous_traj_dir = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data/pocket_trajectory")
candidate_name = [x.stem for x in previous_traj_dir.glob("*")]

if not os.path.exists(affinity_dir):
    os.makedirs(affinity_dir)
# Load the data
with open(pocket_json, "r") as f:
    data = json.load(f)

VE_frcmod = ligand_topology / "VE.frcmod"
VE_lib = ligand_topology / "VE.lib"
mmpbsa_in = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data/VEprotein/mmpbsa.in")

overwrite = True

submit_dir = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data/VEprotein/submit")
if not os.path.exists(submit_dir):
    os.makedirs(submit_dir)

for pocket in data.keys():
    ligand_sdf = result_dir / f"docking.{pocket}.sdf"
    if not os.path.exists(ligand_sdf):
        continue

    protein_name = pocket.split("_")[0]
    if not protein_name in candidate_name:
        continue
    protein_dir = affinity_dir / protein_name
    pocket_dir = protein_dir / pocket
    if os.path.exists(pocket_dir / "FINAL_RESULTS_MMPBSA.dat"):
        continue

    script = [
        "#!/bin/bash",
        "source /etc/profile",
        "module load cuda/11.2",
        "source $HOME/env/amber22.env",
        "set -e",
        "cd /home/gridsan/ywang3/Project/deepmodeling_target_fishing/scripts",
        f"python post_process.py --pocket {pocket}",
    ]
    with open(submit_dir / f"{pocket}.sh", "w") as f:
        f.write("\n".join(script))
    with set_directory(submit_dir):
        run_command(["LLsub", f"{pocket}.sh", "-s", "8"])