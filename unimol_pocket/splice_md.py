from pathlib import Path
import os
import json
import pytraj as pt
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool


def process_traj(case_name):
    print(case_name)
    top = topology_dir/f"{case_name.name}.first_chain.tleap.parm7"
    traj = pt.iterload(str(case_name/"prod.nc"), str(top), frame_slice=[(0, -1, 5)])
    traj.autoimage()
    traj = pt.align(traj, mask="@CA", ref=0)
    out_dir = pocket_dir/case_name.name
    out_dir.mkdir(exist_ok=True)
    # print("write parm7")
    # pt.save(str(out_dir/f"{case_name.name}.parm7"), traj["!:WAT,Na+,Cl-"].top, overwrite=True)
    print("write dcd")
    pt.write_traj(str(out_dir/f"{case_name.name}.dcd"), traj["!:WAT,Na+,Cl-"], overwrite=True)
    os.remove(str(out_dir/f"{case_name.name}.nc"))
    # print("write nc")
    # pt.write_traj(str(out_dir/f"{case_name.name}.nc"), traj["!:WAT,Na+,Cl-"], overwrite=True)
    print("write pdb")
    pt.write_traj(str(out_dir/f"{case_name.name}.pdb"), traj["!:WAT,Na+,Cl-"], frame_indices=[1], overwrite=True)
    return 


data = Path("/home/gridsan/ywang3/Project/deepmodeling_target_fishing/data")
simulation_dir = data / "simulation"
topology_dir = data / "topology"
pocket_dir = data / "pocket_trajectory"


cases = [case for case in simulation_dir.glob("*")]


if __name__ == "__main__":
    with Pool(32) as p:
        p.map(process_traj, cases)

    # get_pocket(case_name)
