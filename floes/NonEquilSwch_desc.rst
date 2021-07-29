* Purpose:

  * This Floe performs Relative Binding Free Energy (RBFE)
    calculations using the Non-Equilibrium Switching (NES) method
    refined by the de Groot lab
    (Gapsys et al., Chem. Sci., 2020, 11, 1140-1152).
* Method Recommendations/Requirements:

  * Four inputs are required:

    * An Orion dataset containing an equilibrium run for each bound ligand
      participating in the NES run.
    * An Orion dataset containing an equilibrium run for each unbound ligand
      participating in the NES run.
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

  * MD, FECalc
* Related Floes:

  * Ligand Bound and Unbound Equilibration for NES [MDPrep] [MD]
  * Equilibrium and Non Equilibrium Switching [MDPrep] [MD] [FECalc]
  * Compare Experimental Affinity with NES Results [Utility] [FECalc]
  * Non-Equilibrium Switching Recovery  [Utility] [FECalc]

The third input above, the text file of edges, describes the map
of desired alchemical transformations of one ligand into another;
each transformation forms an edges of a connected graph of ligands.
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
