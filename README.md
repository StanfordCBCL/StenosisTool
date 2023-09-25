# Probabilistic Intra-Surgical Guidance Tool for Balloon Angioplasty


## Description

We build a pipeline from 3D model to probability analysis, enabling us to answer relevant clinical questions such as optimal path of repair.

## [Paper](NULL)

## Tasks

- [x] Write script for pipeline
- [] Generate Data for Paper
    - [] Regenerate training data
    - [] Rerun Model
    - [] Generate Data to evaluate
- [] Write plotting scripts
    - [] Make sure plots are reproducible
    - [] These can be more hard coded, which is acceptable
- [] Cleanup
    - [] Reorganize code
        - [] Make sure everything is consistent with saving in a config file
        - [] Change svZeroDPlus to pip installed version
        - [] Move svInterface to local instead of submodule
        - [] Comment/Run through to ensure everything works
    - [] Write Documentation
    - [] Create a sample to run
        - [] Ensure all data is present, and config is written

## Installation

The installation is currently only supported for MacOS/Ubuntu.

##### Simvascular

1. Go to [Simvascular](https://github.com/SimVascular/SimVascular) and install either the most recent Github Build, or the official release.
2. Navigate to the Simvascular Home Folder on your local machine and run the shell script `./Simvascular`. You may need to modify the contents depending on your path.

Your Simvascular cmd line should now be set up.

##### Stenosis Tool

1. Clone this repository using `git clone https://github.com/JohnDLee/StenosisTool.git` (Subject to change).
2. If you have conda/miniconda installed, you can create a virtual environment using `conda env install`. Otherwise, install the packages listed in environment.yml
3. You may need to install cuda versions of [Pytorch](https://pytorch.org/get-started/locally/) to run models on Nvidia GPU's
4. Run `pip3 install [-e] .` to install the src library. `[-e]` is for editable mode.

The main stenosis tool should be set up now.

##### svZeroDSolver C++

1. In the Stenosis Tool repository, run `git clone https://github.com/StanfordCBCL/svZeroDPlus.git` (Subject to change).
2. Navigate to the svZeroDPlus directory.
3. Follow instructions under [C++ Documentation](https://stanfordcbcl.github.io/svZeroDPlus/cpp/) to install dependencies and build the C++ Solver

The C++ solver should now be built

apt install libgl1-mesa-glx
sudo apt-get install libxrender1





Write scripts to generate plots, reproducible
Regenerate data.
Retrain model.
Redo data