# OpenMM cubes and workfloes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
* `LigPrepCubes/` - Cubes for preparing molecules
  * `ChargeMCMol` - Assigns partial charges and generates multi-conf molecules with OMEGA
  * `FREDDocking` - Dock MCMols using FRED to a prepared receptor
  * `LigandParameterization` - Parametrize molecules with either GAFF/GAFF2/SMIRNOFF forcefields
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMComplexSetup` - Generate a solvated protein:ligand complex.
  * `OpenMMSimulation` - Runs and OpenMM MD simulation. Minimizes or restarts from saved State.
* `YankCubes/` - YANK cubes
  * `YankHydrationCube` - YANK hydration free energy calculations
  * `YankBindingCube` - YANK absolute binding free energy calculations

## Workfloes
* Testing Floes:
  * `floes/platformTest.py` - Check available OpenMM Platforms
  * `floes/openmm_benchmarking.py` - Performs Benchmarking upon all available Platforms.
* [OpenMM](https://github.com/pandegroup/openmm) Floes:
  * `floes/openmm_MDminimize.py` - Minimize an OpenMM-ready solvated complex
  * `floes/openmm_MDnpt.py` - NPT simulation of an OpenMM-ready solvated complex
  * `floes/openmm_MDnvt.py` - NVT simulation of an OpenMM-ready solvated complex
  * `floes/openmm_MDprep.py` - Set up an OpenMM complex then minimize, warm up and equilibrate a system by using three equilibration stages
  * `floes/openmm_MDprod.py` - Run an unrestrained NPT simulation at 300K and 1atm
* [YANK](https://github.com/choderalab/yank) Floes
  * `floes/yank_hydration.py` - Compute small molecule hydration free energies using YANK.
  * `floes/yank_binding.py` - Compute small molecule absolute binding free energies using YANK.
* Other Floes:
  * `floes/openmm_complex_setup.py` - Setup the protein:ligand complex from **PDBs**.

## Local Installation
```bash
git clone git@github.com:oess/openmm_orion.git
cd openmm_orion

#Create a new local conda environment and install dependencies
conda create -c omnia -c omnia/label/dev -c mobleylab -n dev python=3.5 openmm==7.1.1 openmoltools==0.8.1 ambermini==16.16.0 smarty==0.1.5 parmed==2.7.3 pdbfixer==1.4 smirff99frosst==1.0.5 alchemy==1.2.3 yank==0.15.2
source activate dev

#Install the OpenEye-floe package and toolkits
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
pip install --pre --extra-index-url https://pypi.anaconda.org/OpenEye/channel/beta/simple OpenEye-oenotebook
pip install OpenEye-floe-0.2.158.tar.gz

#Install the main OpenMM Orion Floes
python setup.py develop

# Run the tests.
py.test -v -s PlatformTestCubes
py.test -v -s LigPrepCubes
py.test -v -s -m "not slow" OpenMMCubes
py.test -v -s ComplexPrepCubes
```
