# A Probabilistic Neural Twin for Treatment Planning in Peripheral Pulmonary Artery Stenosis

## Description

The substantial computational cost of high-fidelity models in numerical hemodynamics has, so far, rel-
egated their use mainly to offline treatment planning. New breakthroughs in data-driven architectures
and optimization techniques for fast surrogate modeling provide an exciting opportunity to overcome
these limitations, enabling the use of such technology for time-critical decisions. We discuss an applica-
tion to the repair of multiple stenosis in peripheral pulmonary artery disease through either transcatheter
pulmonary artery rehabilitation or surgery, where it is of interest to achieve desired pressures and flows
at specific locations in the pulmonary artery tree, while minimizing the risk for the patient. Since dif-
ferent degrees of success can be achieved in practice during treatment, we formulate the problem in
probability, and solve it through a sample-based approach. We propose a new offline-online pipeline
for probabilistic real-time treatment planning which combines offline assimilation of boundary condi-
tions, model reduction, and training dataset generation with online estimation of marginal probabilities,
possibly conditioned on the degree of augmentation observed in already repaired lesions. Moreover, we
propose a new approach for the parametrization of arbitrarily shaped vascular repairs through iterative
corrections of a zero-dimensional approximant. 

This Git Directory contains the source code and scripts for replicating our pipeline.

## [Paper (Arxiv)](https://arxiv.org/abs/2312.00854)

## Installation

The installation is currently only supported for MacOS/Ubuntu.

#### Simvascular

Several pre-processing stages of this pipeline was written to connect with Simvascular. The primary requirement is access to `simvascular --python`.

1. Go to [Simvascular](https://github.com/SimVascular/SimVascular) and install either the most recent Github Build, or the official release.
2. Navigate to the Simvascular Home Folder on your local machine and run the shell script `./Simvascular`. You may need to modify the contents depending on your path.

Your Simvascular cmd line should now be set up.

#### Stenosis Tool

1. Clone this repository using `git clone https://github.com/JohnDLee/StenosisTool.git` (Subject to change).
2. If you have conda/miniconda installed, you can create a virtual environment using `conda env create -f environment.yml`. Otherwise, install the packages listed in `environment.yml`
3. You may need to install cuda versions of [Pytorch](https://pytorch.org/get-started/locally/) to run models on Nvidia GPU's. This pipeline was run using `torch==1.13.1+cu117`
4. Run `pip3 install [-e] .` to install the svinterface library. `[-e]` is for editable mode.

#### svZeroDSolver C++

The svZeroDSolver is how 0D simulations are run. There are two options for installation. 

##### Option 1 (pip)
In most cases a pip install is sufficient.

1. `pip install git+https://github.com/simvascular/svZeroDSolver.git`
2. If the building process fails due to missing libraries, they may need to be installed.

##### Option 2 (building using CMake)
The github repository can also by cloned into the main repository and built using CMake.

1. In the Stenosis Tool repository, run `git clone https://github.com/simvascular/svZeroDSolver` (Subject to change).
2. Navigate to the svZeroDPlus directory.
3. Follow instructions under [C++ Documentation](https://simvascular.github.io/svZeroDSolver/) to install dependencies and build the C++ Solver

The setup should now be complete.





    

