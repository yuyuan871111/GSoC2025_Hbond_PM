# Project Management repository for GSoC 2025 - MDAnalysis x ProLIF Project 5: H-Bond interactions from implicit hydrogens

Project outline: [link](https://summerofcode.withgoogle.com/programs/2025/projects/5Otkx8vp)

## Environment installation 
Recommended Python:  3.11 or above
Prerequirement: The python packages are managed with [uv](https://docs.astral.sh/uv/). Please follow the official guide to install `uv`.

Once `uv` is installed, you can clone our repository and synchronize the environment with below code:
```bash
git clone https://github.com/yuyuan871111/GSoC2025_Hbond_PM.git

cd GSoC2025_Hbond_PM

uv sync --python 3.11
```

## Quick start
Look at the `HBond valiadtion - Case study` in this [notebook](./validation/plinder_dataset_validation.ipynb).

For more analyses, you will need to download PLINDER dataset and prepare the protonated structures. See [here](./validation/README.md) for details.
