# Probabilistic Intra-Surgical Guidance Tool for Balloon Angioplasty


## Description

We build a pipeline from 3D model to probability analysis, enabling us to answer relevat clinical questions such as optimal path of repair.

## [Paper](NULL)

## Current Tasks

- [x] Construct Pipeline for Healthy Model

- [] Rewrite Code to be more pipeline-like rather than research
    - Skipped Sections
    - [] Compare 3D 1D, check to make sure its correct.
    - [] Write code to convert 0D to 1D centerlines for visualization
    - [x] Write pre-code for centerline gen
    - [x] Write pre-code for CapInfo????
    - [] Code to compare 3D&0D /  extract 3D to centerlines / 
    - [x] Change BC tuning framework

- [] Construct Pipeline for Unhealthy Model
    - [x] Get Flow File
    - [] Remesh the file so the centerlines work better.
    - Errors
        - [] Code Parasolid -> Polydata does not work in Python. Use GUI
        - [] Code Centerline_Gen does not work in Python. Use GUI.
- [] Write Paper
- [] Clean up extra code.
- [] Create Interface

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