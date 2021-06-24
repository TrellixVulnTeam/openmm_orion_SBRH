* Purpose:

  * This Floe postprocesses a Dataset from protein-ligand MD analysis
    to convert it into a cluster-centric dataset suitable to be used
    as input to other floes.
* Method Recommendations/Requirements:

  * Input is an output dataset from a protein-ligand MD analysis floe:
    * Short Trajectory MD with Analysis
    * Analyze Protein-Ligand MD
* Limitations

  * The resulting cluster-centric dataset is not intrinsically ordered by
    ligand, but this can be done in the spreadsheet of the Analyze page.
* Expertise Level:

  * Regular
* Compute Resource:

  * Minimal
* Keywords:

  * Utility
* Related Floes:

  * Analyze Protein-Ligand MD [MDAnalysis]
  * Short Trajectory MD with Analysis [MDPrep] [MD]

The results datasets output from the protein-ligand MD analysis floes
is ligand-centric, containing cluster information inside.
To be able to select specific clusters to carry forward as input
to another floe, this ligand centric dataset must be converted by this floe
into a cluster-centric dataset.
