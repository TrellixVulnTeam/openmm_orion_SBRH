* Purpose:

  * This Floe performs MD simulations given one or more
    complete solutes as input, each to be treated in its entirety in a separate simulation.
* Method Recommendations/Requirements:

  * The solutes need to have reasonable 3D coordinates, all atoms, and correct
    chemistry (in particular bond orders and formal charges).
  * Each solute can have multiple conformers but each conformer will be run
    separately.
  * The starting configuration should not have very high gradients, in particular
    no bad clashes.
  * Protein components need to be prepared to MD standards: protein chains must be
    capped, all atoms in protein residues (including hydrogens) must be present,
    and missing protein loops resolved or capped.
  * Crystallographic internal waters should be retained where possible.
* Limitations

  * The output produced by this floe is not compatible with the Analyze Protein-Ligand MD floe
  * Currently this floe cannot handle covalent bonds between different components
    such as ligand, protein, and cofactors.
  * Glycosylation on proteins is truncated and the amino acid is capped with H.
* Expertise Level:

  * Regular/Intermediate/Advanced
* Compute Resource:

  * Depends on simulation length; Minimal resources for default 2 ns.
* Keywords:

  * MD, MDPrep
* Related Floes:

  * Bound Protein-Ligand MD [MDPrep] [MD]
  * Ligand Bound and Unbound Equilibration for NES [MDPrep] [MD]
  * Short Trajectory MD with Analysis [MDPrep] [MD]

Each input solute is solvated to make a flask with periodic boundary conditions
and parametrized according to the selected force fields.
A minimization stage is performed on the flask followed
by a warm up (NVT ensemble) and three equilibration stages (NPT ensemble).
In the minimization, warm up and equilibration stages
positional harmonic restraints are applied.
At the end of the equilibration stages a production run (default 2ns)
is performed on the unrestrained flask.
