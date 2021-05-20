The Non Equilibrium Switching (NES) floe performs Relative Binding Affinity Calculations
between a set of provided ligands and the relate provided target protein. The floe
requires also a text file with the map of the edges between the different relative
free energy calculations to run. The file format of the text file is a set of lines
with the syntax:

ligA_name >> ligB_name

where ligA_name and ligB_name are respectively strings of the ligand names for the 
ligand in the starting state A and  the ligand name in the final state B. Because the 
edges to run are defined by ligand names it is important that all the submitted ligands 
have unique names. At the end of the calculations the NES floe will produce a floe 
report where the insight of each edge calculation is reported with different metrics 
used to compute the relative binding affinity.

