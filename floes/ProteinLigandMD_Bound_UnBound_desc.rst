* Purpose:

  * This Floe performs MD simulations given a prepared protein and a set of posed
    and prepared ligands, running both bound and unbound simulations of each ligand
    in preparation for Relative Binding Free Energy Calculations using
    Non-Equilibrium Switching.
* Method Recommendations/Requirements:

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

  * Depends on simulation length; Minimal resources for default 6 ns.
* Keywords:

  * MD, MDPrep
* Related Floes:

  * Bound Protein-Ligand MD [MDPrep] [MD]
  * Short Trajectory MD with Analysis [MDPrep] [MD]

Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately, and the bound and unbound simulations are then carried out.
Currently only one of the conformers will be used in the NES calculations.
A minimization stage is performed on the system followed
by a warm up (NVT ensemble) and several equilibration stages (NPT ensemble).
In the minimization, warm up, and equilibration stages, positional harmonic
restraints are applied on the ligand and protein.
At the end of the equilibration stages a
production run (by default 6 ns) is performed on the unrestrained system.
Two datasets are written, one for the bound and one for the unbound ligands.
