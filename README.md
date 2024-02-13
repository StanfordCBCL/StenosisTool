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

*** Note: Due to the constant updating nature of the svZeroDPlus Solver, the release 1.1 version of svZeroDPlus has been tested to work. The pip installation, `pysvzerod==2.0` is currently being verified ***

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


## Running Pipeline (Overview)

Each tool has been seperated into individual scripts. A more detailed description of what each file does can be found [here](docs/script_documentation.txt). In general, steps are presented assuming a smooth workflow.


### Prerequisites (00)

There are several prerequisites to be satisfied prior to beginning the pipeline.

1. *A prestent VTP (`.vtp`) and corresponding MDL (`.mdl`) file containing the model prior to any stenting procedure. Typically, these can be found in the `Models` directory of a Simvascular project.
2. *N poststent VTP and MDL files containing the model after stenting a SINGLE repair location.
3. An exported CapInfo file from Simvascular, containing the cap name and area for each face in the pre-stent model.
4. A inflow waveform file in the _POSITIVE_ direction, as opposed to the _negative_ direction used in svSolver.
5. (Optionally) RCRT boundary conditions in a (`.rcrt`) file. This file must be in the format accepted by the 0D solver.

*Additionally, the name of the inlet cap for each model must be known (they may be different across models depending on naming convention)


### Setup workspace (01)

We set up an isolated workspace to avoid overwriting information in the original directory.

#### Generating a config file

A config file can be generated using the below command (It is highly recommended to be placed in the same general directory as the data files)

```python3 scripts/01_dev/generate_sv_config.py <outdir>```

The config file looks as follows, with comments on each line corresponding to a prerequisite listed in the prior section:

*** ALL PATHS MUST BE RELATIVE TO THE LOCATION OF THE CONFIG FILE ***

```
# SV Workspace
workspace:
  surface_model: "Models/AS1_SU0308_prestent.vtp"   # pre-stent vtp file (1)
  mdl: "Models/AS1_SU0308_prestent.mdl"             # pre-stent mdl file (1)
  flow_file: "flow_files/inflow_1D.flow"            # inflow waveform file (4)
  rcrt_file:                                        # Optional rcrt file (5) -> if provided, tune parameter should be false
  capinfo: "Models/AS1_SU0308_prestent_CapInfo"     # CapInfo file (3)

# Metadata about model
metadata: 
  model_name: "AS1_SU0308"                          # A name for the model
  diseased: true                                    # diseased flag (typically should be true)
  inlet: "cap_RPA"                                  # the name of the inlet cap to the model

# Options in pipeline
options:
  tune: true                                        # whether the model requires boundary condition tuning -> should be true if rcrt_file is empty

# Tuning Parameters
...                                                 # These can be filled later (and are only necessary if tuning is required).
```

#### Generating workspace

We can generate a workspace anywhere using the command:

```python3 scripts/01_dev/setup_workspace.py -i <path_to_config> -o <path_to_directory>```

Files are copied from their original location and renamed. A new config file is generated, referred to as the main config file from now on.

### Centerline Generation (02)

For each model, pre-stent and post-stent, centerlines must be generated. Depending on the quality of the models, this step may throw an error (centerlines failed to be generated), which can only be resolved by refining the models manually in Simvascular.

Generating the pre-stent centerlines is easily accomplished by running the script below and the centerlines will be tracked automatically.

```python3 scripts/02_centerline_gen/centerline_gen_diseased.py -i <path_to_config>```

Generating the post-stent centerlines for each model must be done manually. It is recommended to construct a directory along the lines of "poststent_centerlines" to save the centerlines. Use

```python3 scripts/02_centerline_gen/centerline_gen_diseased.py -mdl <path_to_poststent_mdl> -vtp <path_to_poststent_vtp> -inlet <name_of_inlet_cap> -o <file_output_path>```

*<file_output_path> is the path to save the output centerlines to and should be of the form `(name).vtp`.

### Lumped Parameter Network (LPN) Setup (03)

#### Generating LPN
With the centerlines, we can use Simvascular's inbuilt functionality to generate a default setup for an LPN/0D model with:

```python3 scripts/03_lpn_setup/lpn_segmentation.py -i <path_to_config>```

This should generate a directory called `LPN_DIR` within the directory containing relevant files. If a `.rcrt` file was provided in the original config, then it is copied into the directory. Otherwise, an empty rcrt with 0'd values is created.

*Note: a workaround was implemented here due to the challenges of connecting Simvascular's Python 3.5 and the newer version of Python, where the actual computation is done in `scripts/03_lpn_setup/sv_lpn_segmentation.py`.

#### Mapping LPN to centerlines
For a future step, it is important to map locations on the LPN and centerlines to each other. Running the script

```python3 scripts/03_lpn_setup/map_junctions_to_centerlines.py```

will map global ids on the centerlines onto the lpn and vice versa.


### Boundary Condition Tuning (if necessary) (04)

The tuning procedure developed here was developed for the limited physiological information provided in our sample dataset and may not be an adequate procedure. Various automated pipelines exist and can replace the one provided here.

The necessary data, which must be filled in the config file, are as follows:

```
tune_params:
  PCWP: 18                  # Pulmonary Capillary Wedge Pressure -> Empirically, this value must be <= to the minPAP value otherwise capacitances may tune to 0.
  R_LPA: 2.81564864366      # *Nonlinear Resistance coefficient for LPA
  R_RPA: 5.5670662746       # *Nonlinear Resistance coefficient for RPA
  mPAP:
  - 42                      # lower bound mean Pulmonary arterial Pressure 
  - 42                      # upper bound mean Pulmonary arterial Pressure 
  maxPAP:
  - 90                      # lower bound systolic/max Pulmonary arterial Pressure 
  - 90                      # upper bound systolic/max Pulmonary arterial Pressure 
  minPAP:
  - 18                      # lower bound diastolic/min Pulmonary arterial Pressure 
  - 18                      # upper bound diastolic/min Pulmonary arterial Pressure 
  rpa_split: 0.52           # The flow split going into the RPA branch
```

*The Nonlinear resistance coefficients can be calculated through the procedure mentioned in Section 3 of the linked paper.

After the data is provided all that needs to be run is:

```python3 scripts/04_tune/tune_bc_nonlinear.py -i <path_to_config>```

### Running 0D Simulations and Visualization

After generating an LPN and determining boundary conditions, solving the LPN is easily accomplished by running:

```python3 scripts/solver_scripts/run_lpn.py -i <config_file> -n <name_for_sim> [-c] [-b] [--l] [--m] [-v]```

```
-n <name_for_sim> # a label for the simulation
-c                # whether to save results as a csv file (recommended)
-b                # whether to convert the c output to the old python output (recommended)
--l               # flag to only output last cycle 
--m               # flag to only output mean values
-v                # whether to output a plot of the last 3 cycles to ensure convergence and verify other features
```

To run the original tuned LPN, the default line below should suffice:

```python3 scripts/solver_scripts/run_lpn.py -i <config_file> -n "no_correction" -c -b --l -v```


### Computing 3D Prestent Ground Truth (05)

To perform optimization of LPNs, we need a ground truth 3D solution to optimize to. Setting up a 3D hemodynamic simulation is built into Simvascular, with instructions found [here](https://simvascular.github.io/documentation/flowsolver.html). The provided scripts convert boundary conditions and compress 3D solutions onto a 1D centerline.

#### Converting 0D RCR to 3D RCR Format
Tuned RCR boundary conditions can be converted to a valid 3D RCR boundary condition file by using the script

```python3 scripts/05_3D_prestent/0D_rcrt_to_3D.py -rcrt <0D_rcrt_file_path> -o <3D_simulation_dir_path>```

where `<0D_rcrt_file_path>` is the path to the rcrt file after tuning (or provided), `<3D_simulation_dir_path>` is the directory containing all the files for a 3D simulation (which MUST include a `.svpre` and `.inp` file).

#### Mapping 3D solution onto 1D Centerlines

`scripts/05_3D_prestent/map_3D_centerlines.py` is a standalone script that can be copied and run on a compute cluster with minimal packages. After uploading the centerlines (which have been mapped to the LPN via `scripts/03_lpn_setup/map_junctions_to_centerlines.py`) to the compute cluster, 3D solutions can be mapped onto those centerlines as such.

```python3 scripts/05_3D_prestent/map_3D_centerlines.py -c <centerlines_vtp_file> -v <3D_solution_volume_vtu_file> -o <output_centerline_vtp_file> <--caps|--juncs|--0D>```

where `<--caps|--juncs|--0D>` are optional mutually exclusive flags to map either all points onto the centerline (>2 hrs) or only the relevant 0D LPN locations (<1 minute).

*** Note: The mapping script has a bug where the flows are computed incorrectly. However, the flows are not used in the current iteration of the pipeline ***

#### Formatting 1D Mapped Solutions 

After mapping solutions onto the 1D centerlines, they should be downloaded and saved to a directory in the workspace. However, the mapping must be formatted one more time via

```python3 scripts/05_3D_prestent/format_3D_centerlines.py -i <config_file> -c <3D_centerline_solution_Vtp> -f <3D_inflow> --s```

where `<3D_inflow>` is the file containing the _negative_ inflow waveform directly used by svSolver to solve the 3D solution.

### Linear Correction


####

####

















