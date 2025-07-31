# GSoC2025_Hbond_PM

Project Management repository for GSoC 2025 - MDAnalysis x ProLIF Project 5:  H-Bond interactions from implicit hydrogens

Project outline: [link](https://summerofcode.withgoogle.com/programs/2025/projects/5Otkx8vp)

## Validation with PLINDER dataset
Implicit Hbond method is currently validated by explicit Hbond method using PLINDER dataset.

### Dataset preparation
First, the receptor (protein) and ligand are protonated using PDB2PQR (PROPKA).
```bash
cd validation

uv run python plinder_dataset_protonate.py
```
Note that several types of errors might occur during the protonation (with PDB2PQR).
1. `cannot convert float NaN to integer`: This happens during hydrogen optimization.
> Current solution: protonate with flag `--noopt`.
* In PLINDER validation set: `1ho4__1__1.B_2.C__1.G_1.H`, `2nws__1__1.A_1.B__1.C_1.D`, `2nws__1__2.A_2.B__2.C_2.D`, `2cth__1__1.A__1.C_1.E_1.F`, `5kod__4__1.D__1.M_1.N`
* In PLINDER test set: `5ayc__1__2.A_3.A__2.B_2.C` (fail), 
`5fux__1__1.A__1.C_1.E_1.F`, `6b63__1__1.A_1.B__1.D_1.E_1.F_1.I_1.L` (fail), `3n5h__1__2.A__2.B_2.D` (fail), `3was__1__2.A_3.A__2.C_2.E` (fail), `7ow2__1__1.A_1.C_1.E_1.G__1.E_1.G`, `2yr6__1__1.A__1.D_1.E_1.F`


2. `Too few atoms present to reconstruct or cap residue`: This happens if there is only 1 heavy atoms in a terminus residue in the topology file.
> Current solution: trim that specific residue and protonate again.
* In PLINDER validation set: `2bo7__2__1.H__1.Y_1.Z`, `3fcs__2__1.C_1.D__1.N`, `5ab0__2__1.B__1.T`, `1o7a__3__2.E__2.K`, `1oi6__1__1.A_1.B__1.D`
* In PLINDER test set: [TODO] `3ek5__1__1.A_1.C_1.D_1.E_1.F__1.K_1.L` (fail), `4og7__1__1.A__1.C`, `2w0w__1__3.A__3.B`, `2iut__1__1.A_1.B__1.E_1.F` (fail), `3bos__1__1.A_1.B__1.M_1.P_1.U` (fail), `1h3c__1__2.B__2.G`, `1h36__1__2.B__2.G`, `1w2x__1__1.A_2.A__2.D` (fail), `2y38__1__1.A__1.B`, `2bs2__1__1.C_1.E_1.F__1.V_1.W` (fail), `2bs3__1__1.C_1.E_1.F__1.V_1.W` (fail)

3. `Biomolecular structure is incomplete:  Found gap in biomolecule structure`: [TODO]
* In PLINDER test set: `6oug__1__1.A_1.B_1.C_1.D__1.I`

4. `?`: [TODO]
* In PLINDER test set: `4c49__3__1.C__1.G`

5. `-9.4 deviates by 0.40000000000000036 from integral, exceeding error tolerance 0.001`: [TODO]
* In PLINDER test set: `1h33__1__1.A__1.D`



### Methodology validation
Then, the comparison and analyses are performed in a jupyter notebook (see details in  `validation/plinder_dataset_validation.ipynb`).