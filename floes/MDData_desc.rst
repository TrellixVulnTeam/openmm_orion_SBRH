* Purpose:

  * This Floe extracts useful results generated from the Analysis floes
    for protein-ligand MD.
* Method Recommendations/Requirements:

  * Required input is a Dataset output by the Analysis floes:

    * Analyze Protein-Ligand MD
    * Short Trajectory MD with Analysis
* Expertise Level:

  * Regular
* Compute Resource:

  * Minimal
* Keywords:

  * Utility
* Related Floes:

  * Analyze Protein-Ligand MD
  * Short Trajectory MD with Analysis

This Floe extracts useful results generated from the Short Trajectory
MD with Analysis floe. The protein, ligand and binding site water trajectories
are extracted with the protein and ligand average and median
clusters. Also the generated HTML Floe report is saved. The MD Data is uploaded
as a .tar.gz file to S3 and it is available under the Files UI tab in Orion
with the user-selected file name.
