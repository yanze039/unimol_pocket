import os
import tqdm
from pathlib import Path
import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from unimol_pocket.logger import getLogger
from unimol_pocket.utils import run_command, set_directory


logger = getLogger()


def grab_first_chain(pdbfile, tag="first_chain", outdir=Path("processed")):
    pdb_name = Path(pdbfile).stem
    u = mda.Universe(pdbfile)
    first_chain_seg_id = u.segments[0].segid
    first_chain = u.select_atoms(f"protein and segid {first_chain_seg_id}")
    if not isinstance(outdir, Path):
        outdir = Path(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    first_chain.write(outdir/f"{pdb_name}.{tag}.pdb")


def PocketPreperation(pdbfile):
    run_command(f"fpocket -f {pdbfile}")


def get_residue_id(pdbfile):

    u = mda.Universe(pdbfile)
    residues = u.select_atoms("protein").residues

    result = []

    for residue in residues:
        residue_id = f"{residue.segid}{residue.resid}"
        result.append(residue_id)

    return result
    

def calc_affinity(mmpbsa_in, cp, rp, lp, crd):
    """Calculate the binding affinity."""
    # MMPBSA.py -O -i mmpbsa.in -cp com.top -rp rec.top -lp lig.top -y traj.crd
    run_command(
        ["MMPBSA.py", "-O", "-i", mmpbsa_in, "-cp", cp, "-rp", rp, "-lp", lp, "-y", crd]
    )
    return


def sdf2pdb(sdf, pdb):
    """Convert sdf to pdb."""
    sdf_stem = Path(sdf).stem
    tmpl = Path(sdf).parent / f"addH_{sdf_stem}.sdf"
    run_command(
        ["obabel", "-i", "sdf", sdf, "-o", "sdf", "-O", tmpl, "-h"]
    )
    run_command(
        ["antechamber", "-i", tmpl, "-fi", "sdf", "-o", pdb, "-fo", "pdb", "-rn", "VE", "-at", "gaff2", "-an", "yes"]
    )



def tleap(script, fname="tleap.in"):
    with open(fname, "w") as fp:
        fp.write(script)
    command = ["tleap", "-f", fname]
    logger.info(f'Running "{" ".join(command)}"')
    return_code, _, _ = run_command(command)
    assert return_code == 0
