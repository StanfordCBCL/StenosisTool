from setuptools import find_packages, setup

# basic setup for pip install -e "."
setup(
    name='svinterface',
    packages=find_packages(),
    install_requires = [
        'scipy >= 1.9.1',
        'numpy >= 1.22.3',
        'vtk >= 9.1.0',
        'pandas >= 1.4.3',
        'matplotlib >= 3.5.0'
        'abc',
        'pyyaml'
    ]
)