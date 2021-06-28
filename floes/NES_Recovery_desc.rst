* Purpose:

  * This Floe allows for the recovery of NES results in case the NES floe hangs.
* Method Recommendations/Requirements:

  * Required input is the Recovery Dataset from the NES floe that hung.
* Expertise Level:

  * Regular
* Compute Resource:

  * Minimal
* Keywords:

  * Utility, FECalc
* Related Floes:

  * Non-Equilibrium Switching [MD] [FECalc]
  * Equilibrium and Non Equilibrium Switching [MDPrep] [MD] [FECalc]

This floe is used to recover the NES results in case the NES floe hangs.
Using as input the Recovery Dataset from NES,
it will finish off the analysis of the NES results to produce
the same output datasets and NES.
The analysis will be based only on the data in the Recovery Dataset.