import os

def collect_pdb(pocket_path, data_path_meta):
    
    """
    collect pdb files from pocket_path to data_path
    pocket_path: output dir of `find_pocket`
    data_path: input dir of `generate_binding_pose`
    """
    
    dirs = os.listdir(pocket_path)
    dirs = [item for item in dirs if item.endswith(".first_chain_out")]
    for dir in dirs:
        # copy .pbd file in dir to ../VEprotein
        outdir = os.path.join(data_path_meta, 'veprotein')
        os.system(f"cp {os.path.join(pocket_path,dir)}/*.pdb {outdir}")
