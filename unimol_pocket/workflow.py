from pathlib import Path
import os
import json
import random
from pathlib import Path
from md import em, heat, pressurize, production
from unimol_pocket.logger import getLogger
from unimol_pocket.utils import set_directory
from unimol_pocket.op import tleap, sdf2pdb, calc_affinity


logger = getLogger()


def maketops(input_dir, output_dir):
    """Create topology files."""

    with set_directory(output_dir):
        for pdb in input_dir.glob("*.pdb"):
            print(pdb)
            try:
                create_simulation_box(
                    protein_pdb=pdb,
                    protein_forcefield="ff19SB",
                    water='tip3p',
                    box_size=15.0,
                    mod_aa=True
                )
            except Exception as e:
                print(e)
                print(f"Failed to create topology for {pdb.stem}")
                raise RuntimeError(f"Failed to create topology for {pdb.stem}")


def run_simulation(
        config,
        pdb_name,
        parm7,
        rst7,
        out_dir
    ):
    """Equilibrate the system."""
    if not isinstance(out_dir, Path):
        out_dir = Path(out_dir)
    pdb_dir = out_dir.joinpath(pdb_name)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if not os.path.exists(pdb_dir):
        os.mkdir(pdb_dir)
    
    tasks = {}
    tasks["path"] = str(pdb_dir.resolve())
    tasks["command"] = []

    with set_directory(pdb_dir):
        logger.info("Pre-heating ...")
        tasks["command"].append(
            heat(
                defname="debug",
                prmtop=parm7, 
                conf=rst7,
                ref=rst7,
                nsteps=5,
                dt=0.00001,
                temp=5.0,
                resstraint_wt=config["heat"]["resstraint_wt"],
                fep=False,
                tempi=4.0,
                ofreq=1,
                fname="debug.in",
                run=False
            )
        )
        tasks["command"].append(
            em(
                "min",
                prmtop=parm7, 
                conf="debug.rst7",
                ref="debug.rst7",
                maxcyc=config["em"]["maxcyc"],
                resstraint_wt=config["em"]["resstraint_wt"],
                fep=False,
                fname="min.in",
                run=False
            ))
        
        logger.info("Pre-heating ...")
        tasks["command"].append(
            
            heat(
                defname="pre_heating",
                prmtop=parm7, 
                conf="min.rst7",
                ref="min.rst7",
                nsteps=20,
                dt=0.0002,
                temp=5.0,
                resstraint_wt=config["heat"]["resstraint_wt"],
                fep=False,
                tempi=4.0,
                ofreq=1,
                fname="pre_heat.in",
                run=False
            ))
        logger.info("Second Emergy minimizing ...")
        tasks["command"].append(
            
            em(
                "min2",
                prmtop=parm7,
                conf="pre_heating.rst7",
                ref="pre_heating.rst7",
                maxcyc=config["em"]["maxcyc"],
                resstraint_wt=config["em"]["resstraint_wt"],
                fep=False,
                fname="min2.in",
                run=False
            ))

        logger.info(f"Heating from 5K to {config['temp']}K ...")
        tasks["command"].append(
            heat(
                defname="heat",
                prmtop=parm7, 
                conf="min2.rst7",
                ref="min2.rst7",
                nsteps=config["heat"]["nsteps"],
                temp=config['temp'],
                resstraint_wt=config["heat"]["resstraint_wt"],
                fep=False,
                ofreq=config["heat"]["ofreq"],
                fname="heat.in",
                run=False
            ))

        logger.info("Pre-Pressurising1 ...")
        tasks["command"].append(
            pressurize(
                defname="pre_press",
                prmtop=parm7, 
                conf="heat.rst7",
                ref="heat.rst7",
                nsteps=3000,
                dt=0.002,
                temp=config['temp'],
                resstraint_wt=config['pressurize_res']['resstraint_wt'],
                irest=1, ntx=5,
                fep=False,
                ofreq=10,
                fname="pre_press.in",
                run=False
            ))

        logger.info("Pre-Pressurising2 ...")
        tasks["command"].append(
            pressurize(
                defname="pre_press2",
                prmtop=parm7, 
                conf="pre_press.rst7",
                ref="pre_press.rst7",
                nsteps=3000,
                dt=0.002,
                temp=config['temp'],
                resstraint_wt=config['pressurize_res']['resstraint_wt'],
                irest=1, ntx=5,
                fep=False,
                ofreq=10,
                fname="pre_press2.in",
                run=False
            ))

        logger.info("Pressurising ...")
        tasks["command"].append(
            pressurize(
                defname="pressurize_res",
                prmtop=parm7, 
                conf="pre_press2.rst7",
                ref="pre_press2.rst7",
                nsteps=config["pressurize_res"]["nsteps"],
                dt=0.002,
                temp=config["temp"],
                resstraint_wt=config["pressurize_res"]["resstraint_wt"],
                irest=1, ntx=5,
                fep=False,
                ofreq=config["pressurize_res"]["ofreq"],
                fname="pressurize_res.in",
                run=False
            ))
        
        logger.info("Pressurising ...")
        tasks["command"].append(
            pressurize(
                defname="pressurize",
                prmtop=parm7, 
                conf="pressurize_res.rst7",
                ref="pressurize_res.rst7",
                nsteps=config["pressurize"]["nsteps"],
                dt=0.002,
                temp=config["temp"],
                resstraint_wt=None,
                irest=1, ntx=5,
                fep=False,
                ofreq=config["pressurize"]["ofreq"],
                fname="pressurize.in",
                run=False
            ))

        logger.info("Production run ...")
        tasks["command"].append( production(
            defname="prod",
            prmtop=parm7, 
            conf="pressurize.rst7",
            ref="pressurize.rst7",
            nsteps=config["production"]["nsteps"],
            dt=0.002,
            temp=config["temp"],
            resstraint_wt=None,
            irest=1, ntx=5,
            iremd=0,
            fep=False,
            ntwe=config["production"]["ntwe"],
            ntwx=config["production"]["ntwx"],
            ntpr=config["production"]["ntpr"],
            ntwr=config["production"]["ntwr"],
            fname="prod.in",
            run=False
        ))
        with open("run.sh", "w") as fp:
            fp.write("\n".join(tasks["command"]))
            fp.write("\n")
            fp.write("touch done.tag")
            fp.write("\n")
    return tasks


def submit_jobs(tasks, env, ngroup=-2, submit=True, n_lambda=-1, out_dir=None):
    if out_dir is None:
        cwd = Path(os.getcwd())
        submit_scripts_dir = cwd.joinpath("submit")
    else:
        submit_scripts_dir = Path(out_dir).joinpath("submit")

    if not os.path.exists(submit_scripts_dir):
        os.mkdir(submit_scripts_dir)

    scripts = {}
    script_head = ["#!/bin/bash"] + env + ["set -e"]

    for task, info in tasks.items():
        script = [f"cd {info['path']}"]
        script += ["bash run.sh"]
        script += ["touch done.tag"]
        scripts[task] = script
    
    if ngroup > 0:
        group_size = len(scripts) // ngroup
        keys = list(scripts.keys())
        random.shuffle(keys)
        for igroup in range(ngroup):
            if igroup < len(scripts) % ngroup:
                igroup_size = group_size + 1
            else:
                igroup_size = group_size

            sub = submit_scripts_dir.joinpath(f"grouped.{igroup}.sh")
            with open(sub, "w") as fp:
                fp.write("\n".join(script_head))
                fp.write("\n")
                for _ in range(igroup_size):
                    # if not os.path.exists(Path(tasks[task]['path']).joinpath("done.tag")):
                    fp.write("\n".join(scripts[keys.pop()]))
                    fp.write("\n")
            if submit:
                with set_directory(submit_scripts_dir):
                    os.system(f"LLsub {sub.name} -s 6 -g volta:1")
        assert len(keys) == 0
    else:
        for task, script in scripts.items():
            sub = submit_scripts_dir.joinpath(f"{task}.sh")
            with open(sub, "w") as fp:
                fp.write("\n".join(script_head))
                fp.write("\n")
                fp.write("\n".join(script))

            if submit:
                # if not os.path.exists(Path(tasks[task]['path']).joinpath("done.tag")):
                with set_directory(submit_scripts_dir):
                    os.system(f"LLsub {sub.name} -s 8 -g volta:1")


MMPBSA_IN = """&general
/
&gb
igb=5, saltcon=0.150,
/

"""


def pocket_scoring(
        ligand_sdf,
        ligand_lib,
        ligand_frcmod,
        protein_pdbfile,
        output_dir
    ):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ligand_name = Path(ligand_sdf).stem
    ligand_pdb = output_dir / f"{ligand_name}.pdb"
    
    sdf2pdb(ligand_sdf, ligand_pdb)
    with open(ligand_pdb, "r") as f:
        txt = f.read()
    with open(ligand_pdb, "w") as f:
        f.write(txt.replace("UNL", "VE "))

    mmpbsa_in = Path(output_dir) / "mmpbsa.in"
    
    with open(mmpbsa_in, "w") as f:
        f.write(MMPBSA_IN)

    with set_directory(output_dir):
        logger.info(f"Create topology files for {ligand_name} ...")
        create_system(protein_pdbfile, "ff14SB", ligand_lib, ligand_frcmod, ligand_pdb)
        logger.info(f"Energy minimization {ligand_name} ...")
        em(
            "min",
            "complex.parm7", 
            "complex.rst7",
            "complex.rst7",
            maxcyc=400,
            fname="min.in",
            run=True,
            cuda=False
        )
        logger.info(f"MMGBSA calculation {ligand_name} ...")
        calc_affinity(mmpbsa_in, cp="complex.parm7", rp="protein.parm7", lp="ligand.parm7", crd="min.rst7")

    
    

def create_system(
                protein_pdb,
                protein_forcefield,
                lib, 
                frcmod, ligand
    ):
    """Create a simulation box with protein and water."""
    if not isinstance(protein_pdb, Path):
        protein_pdb = Path(protein_pdb)

    protein_name = protein_pdb.stem
    protein_pdb = protein_pdb.absolute()

    if os.path.exists(f"{protein_name}.tleap.parm7"):
        logger.info(f"Found existing topology files for {protein_name}. Skip.")
        return
    
    scripts = [
        f"source leaprc.protein.{protein_forcefield}",
        f"source leaprc.protein.{protein_forcefield}_modAA",
        f"loadoff {lib}",
        f"loadamberparams {frcmod}",

        # load the coordinates and create the complex
        f"protein = loadpdb {protein_pdb}",
        # create proteins in solution
        f"ligand = loadpdb {ligand}",
        "complex = combine {protein ligand}",

        # "addions ligands Na+ 0",
        # f"savepdb protein {protein_name}.tleap.pdb",
        f"savepdb complex complex.tleap.pdb",
        f"saveamberparm protein protein.parm7, protein.rst7",
        f"saveamberparm complex complex.parm7, complex.rst7",
        f"saveamberparm ligand ligand.parm7, ligand.rst7",

        "quit"
        
        ]
    fname = f"tleap.{protein_name}.in"
    tleap("\n".join(scripts), fname=fname)




def create_simulation_box(
                        protein_pdb,
                        protein_forcefield,
                        water='tip3p',
                        box_size=15.0,
                        mod_aa=True
    ):
    """Create a simulation box with protein and water."""
    if not isinstance(protein_pdb, Path):
        protein_pdb = Path(protein_pdb)

    protein_name = protein_pdb.stem
    protein_pdb = protein_pdb.absolute()

    if water == "tip3p":
        waterbox = "TIP3PBOX"
        ionparm = "ionsjc_tip3p"
    else:
        logger.error("Water box style can only support tip3p. Not implemented.")
        RuntimeError("Not implemented.")
    
    if os.path.exists(f"{protein_name}.tleap.parm7"):
        logger.info(f"Found existing topology files for {protein_name}. Skip.")
        return
    
    scripts = [
        f"source leaprc.water.{water}",
        f"source leaprc.protein.{protein_forcefield}",
        f"source leaprc.protein.{protein_forcefield}_modAA",
        f"loadAmberParams frcmod.{ionparm}",

        # load the coordinates and create the complex
        f"protein = loadpdb {protein_pdb}",
        # create proteins in solution
        f"solvatebox protein {waterbox} {box_size}",

        # "addions ligands Na+ 0",
        "addionsrand protein Na+ 0",
        f"savepdb protein {protein_name}.tleap.pdb",
        f"saveamberparm protein {protein_name}.tleap.parm7 {protein_name}.tleap.rst7",

        "quit"
        
        ]
    fname = f"tleap.{protein_name}.in"
    tleap("\n".join(scripts), fname=fname)




# if __name__ == "__main__":
#     input_dir = Path(Path("../data/processed").resolve())
#     output_dir = Path(Path("../data/topology").resolve())
#     # maketops(input_dir, output_dir)
    
#     with open("config.json", "r") as fp:
#         config = json.load(fp)
    
#     all_task = {}
#     for pdb_file in output_dir.glob("*.pdb"):
#         pdb_name = pdb_file.stem.split(".")[0]
#         print(pdb_name)
#         parm7 = output_dir.joinpath(f"{pdb_file.stem}.parm7")
#         rst7 = output_dir.joinpath(f"{pdb_file.stem}.rst7")
#         assert os.path.exists(parm7)
#         assert os.path.exists(rst7)
#         tasks = run_simulation(
#             config,
#             pdb_name,
#             parm7,
#             rst7,
#             out_dir=Path(Path("../data/simulation").resolve())
#         )
#         all_task[pdb_name] = tasks
        
#     submit_jobs(all_task, config["env"], ngroup=-2, submit=True, n_lambda=-1, out_dir=Path(Path("../data").resolve()))

#     pass

