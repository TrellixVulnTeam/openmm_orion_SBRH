The *Equilibrium and Non Equilibrium Switching* floe combines equilibrium MD
calculations of the bound and unbound ligands with subsequent
Relative Binding Free Energy calculations using Non-Equilibrium Switching.
Four inputs are required:

* A ligand dataset of prepared ligands posed in the protein active site.
* A protein dataset of the prepared MD-ready protein structure,
  including cofactors and structured waters.
* A text file of edges (explained below), one line per edge,
  of form "ligA_name >> ligB_name"
* [Optional] a text file containing experimental binding free energies
  for at least one ligand, one experimental datapoint per line,
  of form "ligA_name {deltaG(exptl)} {error_deltaG(exptl)} {units}"
  for example, "gn1c -8.56 0.17 kcal/mol".


Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately, and the bound and unbound simulations are then carried out.
The ligands need to have coordinates, all atoms, and correct chemistry. Each
ligand can have multiple conformers but each conformer will be run separately
as a different ligand.
Currently only one of the conformers will be used in the NES calculations.
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

Then, Relative Binding Free Energy (RBFE) calculations are performed
using the Non-Equilibrium Switching (NES) refined by the de Groot lab
(Gapsys et al., Chem. Sci., 2020, 11, 1140-1152).
In addition to the equilibrium trajectories from the previous stage,
additional user input is required here in the form of a
text file of the desired alchemical transformations of
one ligand into another ("edges").
The file must have one line per transformation, of format

ligA_name >> ligB_name

where "ligA_name" and "ligB_name" are the ligand names for
the ligands to be transformed.
These names must correspond exactly to those in
the bound and unbound ligand equilibration datasets.

The floe will draw a number of starting snapshots from
the bound and unbound trajectories,
generate an RBFE alchemical transformation from ligA into ligB,
and carry out the NES fast transformation of ligA into ligB,
and vice versa, for each of the snapshots.
The speed of the NES transformation and the number of snapshots
transformed can be adjusted from default values by the user at runtime.
The floe outputs two floe report/dataset pairs, one for the calculated
RBFE edges (transformations of Ligand A into Ligand B), and
one for the derived predictions of ligand binding free energies.