# Probabilistic Intra-Surgical Guidance Tool for Balloon Angioplasty


## Description

We build a pipeline from 3D model to probability analysis, enabling us to answer relevant clinical questions such as optimal path of repair.

## [Paper](NULL)

## Tasks

- [x] First Pipeline for Healthy Model

- [] Rewrite Code for Organization into pipeline.
    - [] Pipeline
        - [x] Generate config for workspace
        - [x] Setup Workspace
        - [x] Setup LPN 
        - [x] LPN Tuning Setup 
            - [x] !BUG: The LPA and RPA swap??? -> solved by sorting the outlet caps
        - [x] Artificial Stenosis Gen
        - [] Stenosis Detection
        - [] Training Data Generation
        - [] Neural Network Training
        - [] Probability Report
    - [] Research
        - [x] Determine models to use, and set adequate meshing
        - [x] Generate Centerlines for each model (ensure they are correct)
        - [x] Retrieve Cap Info/Inlet information for each model
        - [x] Retrieve Inflow condition for diseased model
        - [] JC vs. Base LPN
            - [] 3D
                - [x] Use 0D boundary conditions to Tune each version
                - [x] Create 3D simulation files for each
                - [X] Run 3D simulation for each
                    - [X] Make it so Diseased model doesn't explode.
                    - [] Make sure Diseased model is solved correctly
                - [] Postprocess results to 1D centerlines to compare
                    - [] Centerline map from 0D
                    - [] Centerline map from 3D
                    - [] Centerline conversion from 3D centerline map to 0D style.
                - [] Ensure values are physiological & values are similar/correlated to 3D
            - [] 0D
                - [] Demonstrate that JC model captures delta P better than Base
                - [] Demonstrate that fixture/AS is insufficient w/ Base model.
    - [] Visualization
        - [x] Better Tuning plots (Add targets)
        - [] Generate representation of 0D model vs 3D model and expansion operation (fixing)
        - [] Interface

- [] Write Paper
- [] Clean Up
    - [] Documentation
    - [] Slides
    - [] Remove extra code.


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