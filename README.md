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
## To get the affinity scores from UniMol predictions
from unimol_pocket.main import pocket_ranking

pocket_ranking(
        ligand_sdf_dir,  # the result dir from Unimol containing the ligand sdf
        output_dir,  # the output path
        protein_dir  # the protein dir containing the processed PDB. 
                     # The same as the pocket output dir.
        )
```