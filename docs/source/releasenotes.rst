#############
Release Notes
#############


v 4.0.1 June 2021
=================

General Notice
---------------------------------------------------------------------------------
* Non-Equilibrium Switching floes for Relative Binding Free Energy (RBFE) affinity calculations
* New knowledge-based ensemble score *BintScore* for pose stability
* Support for taking *Posit* multi-receptor output as input to protein-ligand MD
* OpenFF 1.3.1a1 support
* Bug Fixing

New Floes
--------------------------------------------------------------------------------
* Non-Equilibrium Switching. Relative Binding Free energy  calculations and estimate of Binding Affinity,
  starting from equilibrium MD calculations of the bound and unbound ligands.
* Ligand Bound and Unbound Equilibration for NES. Performs the equilibrium MD calculations of the
  bound and unbound ligands needed for NES, starting from an input protein and a set of posed ligands.
* Equilibrium and Non Equilibrium Switching. Equilibrium MD runs followed by Non-Equilibrium Switching,
  starting from an input protein and a set of posed ligands.
* Non-Equilibrium Switching Recovery. This floe attempts to produce results from uncompleted NES runs,
  starting from the recovery file.
* Compare Experimental Affinity with NES Results. Datasets from a completed NES run can be re-analyzed
  by comparison to a new set of experimental data.

Floe Updates
--------------------------------------------------------------------------------
* Short Trajectory MD with Analysis. The Bint score analysis has been added and *Posit* support
* Analyze Protein-Ligand MD. The Bint score analysis has been added
* Cofactors have been added to the MMPBSA score calculations

New Cubes
--------------------------------------------------------------------------------
* BoundUnboundSwitching cube. This cube accept a flask and emit on the bound port the flask
  containing the a complex while the unbound flask is emitted on the standard output port
* RBFECEdgeGathering cube. This cube gathers the MD equilibrium simulation data related to
  the Bound and Unbound simulations and creates edges used to perform relative binding affinity
  calculations
* NESGMXChimera cube. This cube generates a chimera molecule and injecting it on select frames
  from the Bound and Unbound simulations for the Gromacs forward and reverse Non-Equilibrium simulations
* NESGMX cube.This cubes runs the forward and reverse Bound and Unbound starting equilibrium
  frames by using Gromacs
* NESAnalysis cube. This cube collects the data generated from the Non-Equilibrium runs and calculates
  the relative binding affinities
* PlotNESResults cube. This cube generates plots for the Relative binding affinities and the
  estimated Absolute Affinities versus the provided experimental data. Relevant statistics are calculated
  and tabled
* PredictDGFromDDG cube. This cube applies the Maximum likelihood Estimator method to
  predict the absolute binding affinities starting from the relative binding affinities
* BintScoreInitialPoseAndTrajectory cube. This cube calculates a trajectory ensemble BintScore
  relative to an initial pose.

Cube Updates
--------------------------------------------------------------------------------
* Complex Preparation Cube. This cube has been updated to give better support for clashes
  detected between the ligands and the other flask components for *Posit* support

======================

v3.0.1 November 2020
======================

General Notice
--------------------------------------------------------------------------------
* OpenFF 1.3.0 and 1.2.1 support
* Bug Fixing

======================

v2.0.0 September 2020
======================

General Notice
--------------------------------------------------------------------------------
* Ligand centric cluster Analysis
* New Parametrization Engine
* OpenFF 1.2.0 (Parsley) support
* Bug Fixing

New Floes
--------------------------------------------------------------------------------
* Solvate and Run Protein-Ligand MD. Runs MD on a Protein-Ligand system
* Analyze Protein-Ligand MD. Performs Energetic and Clustering Analysis on the output produced by the
  ``Solvate and Run Protein-Ligand MD`` floe
* Extract Short Trajectory MD results for Download. Extract the results produced by the
  ``Short Trajectory MD with Analysis`` floe as file to download

Floe Updates
--------------------------------------------------------------------------------
* A new Equilibration IV stage has been added to the Short Trajectory MD floe

New Cubes
--------------------------------------------------------------------------------
* A new  ``MD Components`` cube has been added to support the new Parametrization Engine

Cube Updates
--------------------------------------------------------------------------------
* None

======================

v1.0.0
======================

General Notice
--------------------------------------------------------------------------------
* Support for OpenFF 1.1.0 (Parsley), the March 2020 release from the Open Force Field Initiative
* Support for Cuda 10.0
* Exposed a new parameter in the ``TrajToOEMol Cube`` to set the water cutoff
  distance. This distance is used to select water molecules around the protein-ligand
  binding site for each trajectory frame and producing a multi conformer water molecule
* Improved Documentation

New Floes
--------------------------------------------------------------------------------
* None

Floe Updates
--------------------------------------------------------------------------------
* New default parameters have been set on the STMDA floe

New Cubes
--------------------------------------------------------------------------------
* None

Cube Updates
--------------------------------------------------------------------------------
* None

======================


v0.9.6
======================

General Notice
--------------------------------------------------------------------------------
* Support for the two force fields from the Open Force Field Initiative:
  Smirnof99Frosst and OpenFF 1.0.0 (parsley)
* A new automatic color style is applied to the clustering in the
  ``Short trajectory MD analysis`` to match the cluster colors in the Floe Report
* Fixed a wrong setting in the ``Solvation Cube`` that was placing solvent molecules
  too close to the solute. This could have produced un-realistic results for some systems
  where water molecules could have been placed inside proteins
* Fixed a bug in the GAFF/GAFF2 force field where 1-4 interactions were
  not correctly scaled
* Fixed a bug related to un-wanted ligand atom re-ordering
* Fixed a bug in the protein-ligand active site depiction

New Floes
--------------------------------------------------------------------------------
* A new Plain Gromacs floe has been added to perform MD by using Gromacs .tpr files

Floe Updates
--------------------------------------------------------------------------------
* The Yank solvation free energy floe has been removed

New Cubes
--------------------------------------------------------------------------------
* None

Cube Updates
--------------------------------------------------------------------------------
* The force field parametrization cube now support open force field 1.0.0 (parsley)
* The trajectory to multi conformer cube is now adding to the record protein-ligand binding site
  close waters. These are used to perform MMPBSA calculations with explicit water

======================


v0.9.4
======================

General Notice
--------------------------------------------------------------------------------
* Upgrades to the trajectory analysis in floe ``Short Trajectory MD with Analysis`` expose key results to the Analyze page and 3D visualization page in Orion.
* Fixed bug in setting restraints in GROMACS cubes.
* Fixed bug in ordering results in the floe report in floe ``Short Trajectory MD with Analysis``.

New Floes
--------------------------------------------------------------------------------
* The MD Spruce Prep Floe has been removed. Proteins must be prepared with the Spruce Prep floes available in the Classic
   floes now.

* The calculation of MMPBSA can now be also performed by using explicit waters (still experimental)

--------------------------------------------------------------------------------

Floe Updates
--------------------------------------------------------------------------------

* The MD Spruce Prep Floe has been removed

* The Simple MD Floe has been renamed the Plain MD Floe

--------------------------------------------------------------------------------

New Cubes
--------------------------------------------------------------------------------
* A new cube has been developed to check the record size before writing to the Orion backend
    to avoiding floe failures. The new cube has been added to all the floes for sanity check.

* A new cube to estimate the water number around a ligand-protein complex has been developed. The cube is
    used in the MMPBSA calculation with the explicit water flag set on

Cube Updates
--------------------------------------------------------------------------------
* Exposed MMPBSA ensemble average and standard deviation in the :ref:`cube_TrajPBSACube` so that it can be displayed in the Analyze page in Orion.
* :ref:`cube_MDTrajAnalysisClusterReport` now generates trajectory average and median molecules for protein and ligand, with one conformer for each major cluster. These are exposed int the 3D visualization page in Orion.
* :ref:`cube_ClusterOETrajCube` now exposes a link to the per-ligand floe report page so it is available in the Analyze page in Orion.
* In :ref:`cube_ComplexPrepCube` traditional references to the full periodic supermolecular ensemble as a "system" have been replaced with references to a "flask" by analogy with an assay well.
* In :ref:`cube_MDFloeReportCube` the floe report now generates tiled links to individual ligands in the same order as the initial list of ligands.
* In :ref:`cube_MDFloeReportCube` the floe report tiles now show how many major clusters were found for each ligand.
* In :ref:`cube_MDNptCube` and :ref:`cube_MDNvtCube` the restraints are now correctly set in GROMACS for proteins consisting of multiple chains.


* Hint interactions and Styles have been removed from receptors and ligands in the Protein, Ligand and FF parametrization
    setting cubes that could cause problems along the MD analysis stages (debugging is in progress)

* A bug has been fixed in the ligand Elf10 charging cube that was causing problems when carboxylic acid was present
    in a ligand to be charged

* The Trajectory to OEMol, Interaction Energies and PBSA calculation cubes have been updated to account for the explicit
    water in the new MMPBSA calculation

======================
