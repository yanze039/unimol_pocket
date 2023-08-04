from pathlib import Path
import os
from utils import run_command, set_directory
import tqdm

from unimol_pocket.logger import getLogger
from unimol_pocket.op import grab_first_chain, PocketPreperation
from unimol_pocket.workflow import pocket_scoring


logger = getLogger()


MAKE_TOP_SCRIPT = """obabel -i sdf {ligand_name}.sdf  -o sdf -O add_H.sdf -h
antechamber -i add_H.sdf -fi sdf -o mol.gaff2.mol2 -fo mol2 -rn MOL -at gaff2 -an yes -dr no -pf yes -c bcc -nc 0
parmchk2 -i mol.gaff2.mol2 -f mol2 -o mol.frcmod -s gaff2 -a yes
tleap -f - <<_EOF
source leaprc.gaff2
loadamberparams mol.frcmod
VE = loadmol2 mol.gaff2.mol2
saveoff VE mol.lib
quit
_EOF

"""

def find_pocket(pdb_dir, outdir=Path("pocket")):
    if not isinstance(outdir, Path):
        outdir = Path(outdir)
    data_path = Path(pdb_dir)
    pdbs = list(data_path.glob("*.pdb"))
    logger.info("Processing PDB files: remove water and grab the first chain...")
    tag = "first_chain"
    for pdb in tqdm.tqdm(pdbs):
        logger.info(f"Processing {pdb}")
        grab_first_chain(pdb, tag=tag, outdir=outdir)
    
    logger.info("Detecting pockets by Fpocket...")
    for pdb in outdir.glob(f"*.{tag}.pdb"):
        PocketPreperation(pdb)


def pocket_ranking(
        ligand_sdf_dir,
        output_dir,
        protein_dir,
        ligand_name="mol"):
    
    if not isinstance(ligand_sdf_dir, Path):
        ligand_sdf_dir = Path(ligand_sdf_dir)
    
    if not isinstance(output_dir, Path):
        output_dir = Path(output_dir)

    sdfs = ligand_sdf_dir.glob("docking.*.sdf")
    candidate = sdfs[0]
    top_dir = output_dir / "ligand_topology"
    top_dir.mkdir(exist_ok=True)
    
    MAKE_TOP = MAKE_TOP_SCRIPT.format(
        ligand_name=candidate.stem.split(".")[1]
    )
    logger.info("make topology files for the ligand...")
    logger.info(MAKE_TOP)

    with open(top_dir/"maketop.sh", "w") as f:
        f.write(MAKE_TOP)
    
    with set_directory(top_dir):
        run_command(["bash", "maketop.sh"])
    
    lib = top_dir / "mol.lib"
    frcmod = top_dir / "mol.frcmod"
    
    for sdf in sdfs:
        logger.info(f"Scoring sdf file: {str(sdf)}")
        pocket_name = sdf.stem.split(".")[1]
        protein_name = pocket_name.split("_")[0]
        protein_dir = output_dir / protein_name
        protein_dir.mkdir(exist_ok=True)
        pocket_dir = protein_dir / pocket_name
        pocket_dir.mkdir(exist_ok=True)

        pocket_scoring(
                ligand_sdf=sdf,
                ligand_lib=lib,
                ligand_frcmod=frcmod,
                protein_pdbfile=pocket_dir / f"{protein_name}.first_chain.pdb",
                output_dir=pocket_dir
            )
