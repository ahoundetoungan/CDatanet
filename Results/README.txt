Replication Package for Count Data Models with Heterogeneous Peer Effects

This folder contains the necessary files to replicate the results in the main 
paper and the online appendix. 

# R package
I developed the *CDatanet* package, which includes all functions necessary for 
replicating the results. The package can be installed on any operating system 
from CRAN or from GitHub.

# Monte Carlo Simulations
* File `simu-dgpA.R` replicates the simulation results for DGP A (Table 1).
* File `simu-dgpB.R` replicates the simulation results for DGP B (Table 1).
* File `simu-dgpC.R` replicates the simulation results for DGP C (Table 1).
* File  `simu-dgpD.R` replicates the simulation results for DGP D (Table 1).
* File `plot_data.R` plots examples of simulated data for DGPs A, B, C, and D 
(Figure 2).

# Application: Effects of Social Interactions on Participation in Extracurricular 
Activities
* File `0_Add_Health_data.R` prepares the data set to be used and produces the 
data summary in Table D.2 as well as Figure 3.
* File `A_Add_Health_nofixed.R` replicates the estimations without fixed effects 
(Models 1-3 in Tables 2, D.3, and D.4).
* File `B_Add_Health_fixed.R` replicates the estimations without fixed effects 
(Models 4-6 in Tables 2, D.5, and D.6).
* File `C_Add_Health_hete.R` replicates the estimations with heterogeneity in 
peer effects (Models 7 and 8 in Table 2, D.7 and D.8).
* File `D_counterfactual.R` replicates the counterfactual analysis (Figure 4).

## Add Health Data
This research uses data from Add Health, a program directed by Kathleen Mullan 
Harris and designed by J. Richard Udry, Peter S. Bearman, and Kathleen Mullan 
Harris at the University of North Carolina at Chapel Hill, and funded by Grant 
P01-HD31921 from the Eunice Kennedy Shriver National Institute of Child Health 
and Human Development, with cooperative funding from 23 other federal agencies 
and foundations. Information on how to obtain Add Health data files is available 
on the Add Health website (https://addhealth.cpc.unc.edu/). 
