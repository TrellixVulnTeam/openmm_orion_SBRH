* Purpose:

  * This Floe analyzes protein-ligand MD trajectories for pose stability.
* Method Recommendations/Requirements:

  * The input dataset to this floe must be the output dataset from one of the
    other floes that run protein-ligand MD:

    * The "Bound Protein-Ligand MD" floe
    * The "Short Trajectory MD with Analysis" floe (in this case, this floe
      will simply repeat the same analysis).
    * The bound complex dataset from the
      "Ligand Bound and Unbound Equilibration for NES" floe.
* Limitations

  * To avoid excessively large output floe reports, the floe report is
    truncated at the top 100 ligands by ensemble MMPBSA score.
* Expertise Level:

  * Regular
* Compute Resource:

  * Minimal
* Keywords:

  * MDAnalysis
* Related Floes:

  * Short Trajectory MD with Analysis [MDPrep] [MD] [MDAnalysis]

Trajectories from different starting poses of the same ligand are combined and
analysed collectively.
One analysis is in terms of interactions between the
ligand and the active site.
Another looks at clustering the ligand positions in the protein active site
after fitting the trajectory based on active site C_alphas.
Ensemble MMPBSA and ensemble BintScore calculations are carried out
on the trajectory and are localized to the ligand clusters.
An HTML Floe report is generated for the top 100 ligands by ensemble MMPBSA score.
