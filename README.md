Summer Project for John Lee

- Stenosis Tool

To set up project from root

Install a version of Simvascular from June 2022

For python requirements:

    1. conda env create --name <env_name> --file environment.yml
    2. pip install git+https://github.com/JohnDLee/svZeroDSolver.git (Slightly modified version of original svZeroDSolver)
    3. pip install -e .

For cpp svZeroDSolver:

    1. git clone https://github.com/richterjakob/svZeroDSolver.git
    2. cd svZeroDSolver
    3. If cmake is not installed, run "brew install cmake"
    4. Follow instructions at https://richterjakob.github.io/svZeroDSolver/cpp/ to install dependencies & cmake

