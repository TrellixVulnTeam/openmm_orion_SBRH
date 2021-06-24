* Purpose:

  * This Floe performs short MD simulations given a prepared protein and a set of posed
    and prepared ligands, then analyzes the trajectory for pose stability.
* Method Recommendations/Requirements:

  * To avoid excessively large output floe reports, the floe report is
    truncated at the top 100 ligands by ensemble MMPBSA score.
  * The ligands need to have reasonable 3D coordinates, all atoms, and correct
    chemistry (in particular bond orders and formal charges).
  * Each ligand can have multiple conformers but each conformer will be run
    separately as a different ligand.
  * The starting poses should not have very high gradients, in particular
    no bad clashes with the protein.
  * The protein needs to be prepared to MD standards: protein chains must be
    capped, all atoms in protein residues (including hydrogens) must be present,
    and missing protein loops resolved or capped.
  * Crystallographic internal waters should be retained where possible.
* Limitations

  * Currently this floe cannot handle covalent bonds between different components
    such as ligand, protein, and cofactors.
  * Glycosylation on proteins is truncated and the amino acid is capped with H.
* Expertise Level:

  * Regular/Intermediate/Advanced
* Compute Resource:

  * Minimal
* Keywords:

  * MD, MDPrep, MDAnalysis
* Related Floes:

  * Bound Protein-Ligand MD [MDPrep] [MD]

    * First half of this floe
  * Analyze Protein-Ligand MD [MDAnalysis]

    * Last half of this floe
  * Convert MD Analysis results to Cluster-Centric Dataset [Utility]

    * Convert ligand-centric output from this floe into cluster-centric
      output to select clusters for further work
  * Extract Short Trajectory MD Results for Download [Utility]

    * Extract and save in a .tar.gz file:

      * The protein, ligand and binding site water trajectories as
        multi-conformer OEMols.
      * The Average and Median protein-ligand complex for each cluster.

Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately, and the complex is solvated and parametrized according to the
selected force fields. A minimization stage is performed on the system followed
by a warm up (NVT ensemble) and three equilibration stages (NPT ensemble). In the
minimization, warm up, and equilibration stages, positional harmonic restraints are
applied on the ligand and protein. At the end of the equilibration stages a short
(default 2ns) production run is performed on the unrestrained system.
The production run is then analyzed.
Trajectories from different starting poses of the same ligand are combined and
analysed collectively.
One analysis is in terms of interactions between the
ligand and the active site.
Another looks at clustering the ligand positions in the protein active site
after fitting the trajectory based on active site C_alphas.
Ensemble MMPBSA and ensemble BintScore calculations are carried out
on the trajectory and are localized to the ligand clusters.
An HTML Floe report is generated for the top 100 ligands by ensemble MMPBSA score.
