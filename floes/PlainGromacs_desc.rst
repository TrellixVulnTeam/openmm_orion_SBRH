* Purpose:

  * This Floe performs MD simulations using *GROMACS* given
    ready-to-go input *.tpr* files.
* Method Recommendations/Requirements:

  * The input *.tpr* files need to be correctly prepared (no checking is done).
* Limitations

  * The output produced by this floe is not compatible with the Analyze Protein-Ligand MD floe.
* Expertise Level:

  * Advanced: this floe is for *GROMACS* MD experts.
* Compute Resource:

  * Depends on simulation length as given in the *.tpr* input file.
* Keywords:

  * MD
* Related Floes:

  * Solvate and Run MD [MDPrep] [MD]
  * Bound Protein-Ligand MD [MDPrep] [MD]
  * Ligand Bound and Unbound Equilibration for NES [MDPrep] [MD]

The Floe will run the MD in stages, where each stage runs for *n* hours
(10hrs default) in between checkpointing in a recovery dataset.
The trajectory accumulated in each stage is stored in gromacs file format
in the user's filespace.
The Gromacs Cube then
restarts from the recovery dataset and runs for
an additional *n* hours in a cycle
until completing the number of md steps specified in
the input *.tpr* file.
If the recovery dataset is
provided as input, Gromacs will start from the last
checkpoint saved in the recovery dataset. If both *.tpr* 
and recovery dataset files are provided the recovery dataset 
will supersede the *.tpr* file.
