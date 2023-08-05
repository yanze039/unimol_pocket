from pathlib import Path
import MDAnalysis as mda
import os
import json

def get_pocket_json(pocket_path, data_path_meta):
    
    """
    generate json file of residue id of each pocket for each protein
    whose keys are protein_pocketid, and values are list of residue id

    pocket_path: output dir of `find_pocket`
    data_path_meta: input dir of `generate_binding_pose`
    """
    
    protein_collect = os.listdir(pocket_path)
    protein_collect = list(set([item[:4] for item in protein_collect]))
    pocketdict = {}
    for protein in protein_collect:
        pocketslist = os.listdir(os.path.join(pocket_path, f"{protein}.first_chain_out/pockets"))
        pocketslist = [item for item in pocketslist if item.endswith(".pdb")]
        for pocket in pocketslist:
            pid = pocket[-9]
            pdb_file = Path(os.path.join(pocket_path, f"{protein}.first_chain_out/pockets/{pocket}"))
            u = mda.Universe(pdb_file)
            residues = u.select_atoms("protein").residues
            result = []
            for residue in residues:
                residue_id = f"{residue.segid}{residue.resid}"
                result.append(residue_id)
            result = list(set(result))
            pocketdict[f'{protein}_{pid}'] = result

    json_object = json.dumps(pocketdict, indent=4)
    with open(os.path.join(data_path_meta, "veprotein.pocket.json"), "w") as outfile:
        outfile.write(json_object)
