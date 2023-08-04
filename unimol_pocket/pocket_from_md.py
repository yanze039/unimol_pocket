from pathlib import Path
import os
from multiprocessing import Pool
from contextlib import contextmanager
from pathlib import Path
import subprocess
from typing import Optional, List


# def process_traj(case_name):
#     print(case_name)
#     top = topology_dir/f"{case_name.name}.first_chain.tleap.parm7"
#     traj = pt.iterload(str(case_name/"prod.nc"), str(top))
#     traj.autoimage()
#     traj = pt.align(traj, mask="@CA", ref=0)
#     out_dir = pocket_dir/case_name.name
#     out_dir.mkdir(exist_ok=True)
#     # pt.write_traj(str(lmd_dir.joinpath("prod.new.pdb")), traj["!:WAT,Na+,Cl-"], options='model', frame_indices=[i*100 for i in range(10)], overwrite=True)
#     # print("write parm7")
#     # pt.save(str(out_dir/f"{case_name.name}.parm7"), traj["!:WAT,Na+,Cl-"].top, overwrite=True)
#     # # print("write dcd")
#     # # pt.write_traj(str(out_dir/f"{case.name}.dcd"), traj["!:WAT,Na+,Cl-"], overwrite=True)
#     # print("write nc")
#     # pt.write_traj(str(out_dir/f"{case_name.name}.nc"), traj["!:WAT,Na+,Cl-"], overwrite=True)
#     print("write pdb")
#     pt.write_traj(str(out_dir/f"{case_name.name}.pdb"), traj["!:WAT,Na+,Cl-"], frame_indices=[1], overwrite=True)
#     return 

@contextmanager
def set_directory(path: Path):
    """Sets the current working path within the context.
    Parameters
    ----------
    path : Path
        The path to the cwd
    Yields
    ------
    None
    
    Examples
    --------
    >>> with set_directory("some_path"):
    ...    do_something()
    """
    cwd = Path().absolute()
    path.mkdir(exist_ok=True, parents=True)
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)


def run_command(
        cmd: List, 
        stdin: Optional[str] = None,
        shell: Optional[bool] = None,
        print_to_terminal: bool = True
    ):
    
    if print_to_terminal:
        stdout = None
    else:
        stdout = subprocess.PIPE
    
    stderr = subprocess.STDOUT
    
    with subprocess.Popen(
        args=cmd,
        stdin=subprocess.PIPE,
        stdout=stdout,
        stderr=stderr,
        encoding="utf-8",
        shell=shell
    ) as subp:
        out, err = subp.communicate(input=stdin)
        return_code = subp.poll()
    print(out)
    print(err)
    return return_code, out, err

data = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data")
simulation_dir = data / "simulation"
topology_dir = data / "topology"
pocket_dir = data / "pocket_trajectory"


cases = [case for case in simulation_dir.glob("*")]


def get_pocket(case_name):
    # mdpocket --trajectory_file trajectory_superimposed.dcd --trajectory_format dcd -f reference.pdb
    top = topology_dir/f"{case_name.name}.first_chain.tleap.parm7"
    out_dir = pocket_dir/case_name.name
    ref = (out_dir / f"{case_name.name}.pdb")
    traj = (out_dir / f"{case_name.name}.dcd")
    fm = "dcd"
    # os.chdir()
    print(f"mdpocket --trajectory_file {traj.name} --trajectory_format {fm} -f {ref.name}")
    print(out_dir)
    with set_directory(out_dir):
        run_command(f"mdpocket --trajectory_file {traj.name} --trajectory_format {fm} -f {ref.name}".split())
    # run_command(f"mdpocket --trajectory_file {str(traj)} --trajectory_format {fm} -f {str(ref)}".split())


if __name__ == "__main__":
    with Pool(32) as p:
        p.map(get_pocket, cases)
    
