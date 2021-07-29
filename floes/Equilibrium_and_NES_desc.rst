* Purpose:

  * This Floe performs Relative Binding Free Energy (RBFE)
    calculations using the Non-Equilibrium Switching (NES) method
    refined by the de Groot lab
    (Gapsys et al., Chem. Sci., 2020, 11, 1140-1152).
    It also carries out the equilibrium MD runs which must precede NES.
* Method Recommendations/Requirements:

  * Four inputs are required:

    * A protein prepared to MD standards: protein chains must be
      capped, all atoms in protein residues (including hydrogens) must be present,
      and missing protein loops resolved or capped.
      Crystallographic internal waters should be retained where possible.
    * A dataset of posed ligands need to have reasonable 3D coordinates,
      all atoms, and correct
      chemistry (in particular bond orders and formal charges).
      The starting poses should not have very high gradients, in particular
      no bad clashes with the protein.
    * A text file of edges (explained below), one line per edge,
      of form "ligA_name >> ligB_name".
    * [Optional] a text file containing experimental binding free energies
      for at least one ligand, one experimental datapoint per line,
      of form "ligA_name {deltaG(exptl)} {error_deltaG(exptl)} {units}"
      for example, "gn1c -8.56 0.17 kcal/mol".
* Limitations

  * If no experimental binding free energies (the fourth input above)
    are given, the estimation of ligand binding free energies has no
    reference value so the relative values will be centered
    around the mean.
  * Currently there is no mitigation for the effects of changes in
    buried waters, protein sidechain flips, or large protein movements
    between ligA and ligB.
* Expertise Level:

  * Regular
* Compute Resource:

  * High
* Keywords:

  * MDPrep, MD, FECalc
* Related Floes:

  * Ligand Bound and Unbound Equilibration for NES [MDPrep] [MD]
  * Non Equilibrium Switching [MD] [FECalc]
  * Compare Experimental Affinity with NES Results [Utility] [FECalc]
  * Non-Equilibrium Switching Recovery  [Utility] [FECalc]

This floe combines equilibrium MD
calculations of the bound and unbound ligands with subsequent
Relative Binding Free Energy calculations using Non-Equilibrium Switching.
Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately, and the bound and unbound simulations are then carried out.
Each ligand can have multiple conformers but each conformer will be run
separately as a different ligand.
Currently only one of the conformers will be used in the NES calculations.
A minimization stage is performed on the system followed
by a warm up (NVT ensemble) and several equilibration stages (NPT ensemble).
In the minimization, warm up, and equilibration stages, positional harmonic
restraints are applied on the ligand and protein.
At the end of the equilibration stages a
production run (by default 6 ns) is performed on the unrestrained system.
Separate datasets are written for the bound and unbound ligands.

Then, Relative Binding Free Energy (RBFE) calculations are performed
using the Non-Equilibrium Switching (NES) method.
Here the third input mentioned above is used, a
the text file of edges, describing the map
of desired alchemical transformations of one ligand into another.
Each transformation forms an edges of a connected graph of ligands.
The file must have one line per transformation, of format

ligA_name >> ligB_name

where "ligA_name" and "ligB_name" are the ligand names for
the ligands to be transformed.
These ligand names must correspond exactly to those in
the bound and unbound ligand equilibration datasets.

The floe will draw a number of starting snapshots from
the bound and unbound trajectories of the ligands.
Then for each edge in the edge file, it will
generate an RBFE alchemical transformation from ligA into ligB,
and carry out the NES fast transformation of ligA into ligB,
and vice versa, for each of the snapshots.
The resulting relative free energy change, or DeltaDeltaG,
for each edge is the primary output of this method.
A maximum likelihood estimator is then used to derive
a predicted binding affinity (free energy, or DeltaG) for each ligand.
The mean value of the input experimental binding free energies
is used as the reference value for the computed ones.

The speed of the NES transformation and the number of snapshots
transformed can be adjusted from default values by the user at runtime.
The floe outputs two floe report/dataset pairs, one for the calculated
RBFE edges (DeltaDeltaGs), and
one for the derived affinity predictions (DeltaGs) of ligand.
