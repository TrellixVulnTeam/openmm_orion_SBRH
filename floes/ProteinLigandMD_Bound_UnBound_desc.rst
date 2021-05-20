This floe carries out the equilibrium MD simulations of the ligand, bound and
unbound, for Non-Equilibrium Switching relative binding free energy calculations.
Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately, and the bound and unbound simulations are then carried out.
The ligands need to have coordinates, all atoms, and correct chemistry. Each
ligand can have multiple conformers but each conformer will be run separately
as a different ligand.
The protein needs to be prepared to MD standards: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved. Crystallographic internal waters should be retained where
possible. The parametrization of some common nonstandard residues is partially
supported.
A minimization stage is performed on the system followed
by a warm up (NVT ensemble) and several equilibration stages (NPT ensemble).
In the minimization, warm up, and equilibration stages, positional harmonic
restraints are applied on the ligand and protein.
At the end of the equilibration stages a
production run (by default 6 ns) is performed on the unrestrained system.
Separate datasets are written for the bound and unbound ligands.

Required Input Parameters:
* A ligand dataset of prepared ligands posed in the protein active site.
* A protein dataset of the prepared MD-ready protein structure, including cofactors and structured waters.
