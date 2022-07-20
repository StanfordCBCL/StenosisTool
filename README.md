Summer Project for John Lee

- Stenosis Tool

To set up project from root

For python requirements:

    conda env create --name <env_name> --file environment.yml
    pip install git+https://github.com/JohnDLee/svZeroDSolver.git (Slightly modified version of original svZeroDSolver)
    pip install -e .[dev] (".[dev]" for Macos)


For cpp svZeroDSolver:

    1. git clone https://github.com/richterjakob/svZeroDSolver.git
    2. cd svZeroDSolver
    3. If cmake is not installed, run "brew install cmake"
    4. Follow instructions at https://richterjakob.github.io/svZeroDSolver/cpp/ to install dependencies & cmake

