# A pocket ranking pepline based on UniMol

## Environment Setup

* FPocket is required to detect the pockets. For installation instruction, please refer to https://github.com/Discngine/fpocket.

* Amber/AmberTools is required to conduct the MMPBSA calculations. For installation instruction, please refer to https://ambermd.org/Installation.php.

* UniMol is required to predict the binding pose. For installation instruction, please refer to https://github.com/dptech-corp/Uni-Mol.

## Usage Example

```python
## To find the pockets for PDBs in a folder
from pathlib import Path
from unimol_pocket.main import find_pocket

find_pocket(pdb_dir, # the protein dir containing raw PDB
            outdir=Path("pocket"))
# A output folder will be generated.
```

```python
## To generate binding pose with UniMol
from unimol_pocket.main import binding_pose

binding_pose(
        data_path_meta, # data path of binding pose prediction input files
        weight_path, # dir of the pretrained model
        smiles_list=['CC1=C2C(OC(CC)(C)CC2)=C(C)C(C)=C1O'], # drugs of interest, VE as an example
        pocket_path=Path("pocket") # output dir of `find_pocket`
        )
# A output folder containing ligand sdf will be generated.
```

```python
## To get the affinity scores from UniMol predictions
from unimol_pocket.main import pocket_ranking

pocket_ranking(
        ligand_sdf_dir,  # the result dir from Unimol containing the ligand sdf
        output_dir,  # the output path
        protein_dir  # the protein dir containing the processed PDB. 
                     # The same as the pocket output dir.
        )
```
