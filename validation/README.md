## Validation with PLINDER dataset
Implicit H-bond method is currently validated by explicit H-bond method using [PLINDER](https://www.plinder.sh/)  validation (n ~= 1100) and testing (n ~= 1400) datasets (v2).

Here is the workflow to reproduce our validation:
1. [Dataset preparation](#1-dataset-preparation).
    > If you don't want to prepare the files yourself, we have data archive for the PLINDER validation and testing set [here](https://doi.org/10.5281/zenodo.16950708) or in [mirror archive](https://drive.google.com/drive/folders/1ZZ9zrzBI82vunOM4xmkiUtQX5Sm-eqND) (`plinder_protonated_systems.zip`).
2. [Calculating H-bond interactions](#2-calculating-h-bond-interactions).
    > We have precomputed data [here](https://doi.org/10.5281/zenodo.16950708) or in [mirror archive](https://drive.google.com/drive/folders/1ZZ9zrzBI82vunOM4xmkiUtQX5Sm-eqND) (Filenames: `plinder_val_fps.zip` and `plinder_test_fps.zip`. It is noted that some systems are unavailable due to the failure of reading the topology file.).
3. [Setting the thresholds and testing](#3b-setting-the-thresholds-for-the-implicit-h-bond-methods-using-plinder-validation-set-test-the-performance-on-plinder-testing-set).
    > If you are using the precomputed data, make sure the data are in the correct file path.

### 1. Dataset preparation
_You can skip this step if you use our data archive._
First of all, the receptor (here, protein) and ligand will need to be downloaded to your local machine. Learn how to download the data from [PLINDER documentation](https://www.plinder.sh/).

Then, receptors are protonated using `PDB2PQR` (`PROPKA`), and hydrogens are added to ligands using openbabel. You can simply run the below scripts:
```bash
cd validation

uv run python plinder_dataset_protonate.py # for the PLINDER validation/test set
# will take ~3 hrs using one thread
```
Note that if you want to protonate your structures in the validation set, please modify the path of dataframe in `plinder_dataset_protonate.py`. For example:
```python
...

df = pd.read_csv("./data/plinder_val_systems.csv") # PLINDER validation set

...
```
It is noted that the dataset preparation might face a couple of issues, we dealt with those issues by either slightly modifying topology or replacing `PDB2PQR` with `reduce` (another software to add hydrogens) (see details in [Troubleshooting - Dataset preparation](#plinderpinder-dataset-preparation-protonation)).

Once the dataset preparation finished, you can follow our steps to compute the H-bond interactions using implicit and explicit methods.

### 2. Calculating H-bond interactions
_You can skip this step if you use our data archive._
```bash
cd validation

mkdir -p ./val/explicit
mkdir -p ./val/implicit
mkdir -p ./test/explicit
mkdir -p ./test/implicit

uv run python plinder_dataset_fp_calculation.py # for the PLINDER validation set
# will take ~3 hrs using one thread
uv run python plinder_dataset_fp_calculation_for_testing.py # for the PLINDER testing set
# will take ~3 hrs using one thread
```
### 3a. One example from the PLINDER dataset for calculating the implciit and explicit H-bond interactions 
* Please looks at [this notebook](./plinder_dataset_validation_single.ipynb). This notebook is 80% simlar to the [simple example](../implicit_vs_explicit.ipynb). The only difference is that we demonstrate how to calculate H-bond interactions if the protein interacting with multiple ligands (basically, iteratively load the ligand and calculate fingerprints in sequence).  

### 3b. Setting the thresholds for the implicit H-bond methods using PLINDER validation set (test the performance on PLINDER testing set)
_If you use our data archive, please make sure the file path is correct._
* To decide the thresholds (for `acceptor/donor atom angles' deviation`, `acceptor/donor plane angles`, and `autodock_vina_hbond_potential` rescaling factors), we performed the comparison and analyses in [this notebook](./plinder_dataset_validation.ipynb).
* To test the performance when applying new thresholds, please see [this notebook](./plinder_dataset_test.ipynb).


## Validation with PINDER dataset (protein-protein interaction)
Similar to the validation with PLINDER dataset, we followed a similar protocol using [PINDER](https://www.pinder.sh/) dataset (2024-02).

### 1. Dataset preparation
Considering it is a protein-protein dataset, it is no need to use openbabel to add hydrogens. We used `PDB2PQR` (`PROPKA`) and `reduce` to prepare the structures.

If you don't want to prepare the files yourself, we have data archive for the PINDER validation and testing set [here](https://doi.org/10.5281/zenodo.16950708) or in [mirror archive](https://drive.google.com/drive/folders/1ZZ9zrzBI82vunOM4xmkiUtQX5Sm-eqND) (`pinder_protonated_pdbs.zip`). 
```bash
cd validation

uv run python pinder_dataset_protonate.py # for the PINDER validation set
# will take ~3 hrs using one thread
uv run python pinder_dataset_protonate_for_testing.py # for the PINDER testing set
# will take ~3 hrs using one thread
```


### 2. Calculating H-bond interactions
You can skip this step if using our precomputed data, which can be found [here](https://doi.org/10.5281/zenodo.16950708) or in [mirror archive](https://drive.google.com/drive/folders/1ZZ9zrzBI82vunOM4xmkiUtQX5Sm-eqND) (Filenames: `pinder_val_fps.zip` and `pinder_test_fps.zip`. It is noted that some systems are unavailable due to the failure of reading the topology file.).
```bash
cd validation

# create the output folder
mkdir -p ./pinder_val/explicit
mkdir -p ./pinder_val/implicit
mkdir -p ./pinder_test/explicit
mkdir -p ./pinder_test/implicit

uv run python pinder_dataset_fp_calculation.py # for the PINDER validation set
# will take ~3 hrs using one thread
uv run python pinder_dataset_fp_calculation_for_testing.py # for the PINDER testing set
# will take ~3 hrs using one thread
```

### 3. Setting the thresholds for the implicit H-bond methods using PINDER validation set (test the performance on PINDER testing set)
_If you are using the precomputed data, make sure the data are in the correct file path._
* To decide the thresholds, we performed the simliar analyeses using PINDER validation set in [this notebook](./pinder_dataset_validation.ipynb).
* To test the performance when applying new thresholds, please look at [this notebook](./pinder_dataset_test.ipynb).



## Troubleshooting
### PLINDER/PINDER Dataset preparation (protonation)
Note that several types of errors might occur during the protonation (with PDB2PQR).
#### 1. `cannot convert float NaN to integer`: This happens during hydrogen optimization.
> Current solution: Try to protonate with flag `--noopt`. If still not working, use "reduce" instead (however, the water molecules will be removed and HIS is always HIP).
* In PLINDER validation set: `1ho4__1__1.B_2.C__1.G_1.H`, `2nws__1__1.A_1.B__1.C_1.D`, `2nws__1__2.A_2.B__2.C_2.D`, `2cth__1__1.A__1.C_1.E_1.F`, `5kod__4__1.D__1.M_1.N`
* In PLINDER test set: `5ayc__1__2.A_3.A__2.B_2.C` (with reduce), 
`5fux__1__1.A__1.C_1.E_1.F`, `6b63__1__1.A_1.B__1.D_1.E_1.F_1.I_1.L` (with reduce), `3n5h__1__2.A__2.B_2.D` (with reduce), `3was__1__2.A_3.A__2.C_2.E` (with reduce), `7ow2__1__1.A_1.C_1.E_1.G__1.E_1.G`, `2yr6__1__1.A__1.D_1.E_1.F`.
* IN PINDER validation set: ... (some)


#### 2. `Too few atoms present to reconstruct or cap residue`: This happens if there is only 1 heavy atoms in a terminus residue in the topology file.
> Current solution: Try to trim that specific residue and protonate again. If still not working, use "reduce" instead (however, the water molecules will be removed and HIS is always HIP).
* In PLINDER validation set: `2bo7__2__1.H__1.Y_1.Z`, `3fcs__2__1.C_1.D__1.N`, `5ab0__2__1.B__1.T`, `1o7a__3__2.E__2.K`, `1oi6__1__1.A_1.B__1.D`
* In PLINDER test set: `3ek5__1__1.A_1.C_1.D_1.E_1.F__1.K_1.L` (with reduce), `4og7__1__1.A__1.C`, `2w0w__1__3.A__3.B`, `2iut__1__1.A_1.B__1.E_1.F` (with reduce), `3bos__1__1.A_1.B__1.M_1.P_1.U` (with reduce), `1h3c__1__2.B__2.G`, `1h36__1__2.B__2.G`, `1w2x__1__1.A_2.A__2.D` (with reduce), `2y38__1__1.A__1.B`, `2bs2__1__1.C_1.E_1.F__1.V_1.W` (with reduce), `2bs3__1__1.C_1.E_1.F__1.V_1.W` (with reduce).
* IN PINDER validation set: `2jbx__A1_Q85295-R.pdb`, `2jbx__B1_Q85295-L.pdb`, `2wmp__A1_Q1R2W9-R.pdb`, ... (a lot)

#### 3. `Biomolecular structure is incomplete:  Found gap in biomolecule structure`:  This happens during biomolecule debumping.
> Current solution: Try to protonate with flag `--noopt`. If still not working, use "reduce" instead (however, the water molecules will be removed and HIS is always HIP).
* In PLINDER test set: `6oug__1__1.A_1.B_1.C_1.D__1.I` (with reduce)
* IN PINDER validation set: ... (some)

#### 4. `ZeroDivisionError: float division by zero in PROPKA`: This happens during PROPKA calculation.
> Current solution: Use "reduce" instead (however, the water molecules will be removed and HIS is always HIP).
* In PLINDER test set: `4c49__3__1.C__1.G` (with reduce).
* IN PINDER validation set: ... (some)

#### 5. `-9.4 deviates by 0.40000000000000036 from integral, exceeding error tolerance 0.001`: This happens when applying force field to biomolecule states.
> Current solution: Try to protonate with flag `--noopt`. If still not working, use "reduce" instead (however, the water molecules will be removed and HIS is always HIP).
* In PLINDER test set: `1h33__1__1.A__1.D` (with reduce).
* IN PINDER validation set: ... (some)