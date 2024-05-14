**Replication codes for *[Count Data Models with Heterogeneous Peer Effects under Rational Expectations](https://dx.doi.org/10.2139/ssrn.3721250)***

### Simulations
- `simu-dgpA.R` replicates the simulation results for DGP A.
- `simu-dgpB.R` replicates the simulation results for DGP B.
- `simu-dgpC.R` replicates the simulation results for DGP C.
- `simu-dgpD.R` replicates the simulation results for DGP D.
- `plot_data.R` plots examples of simulated data for DGPs A, B, C, and D (Figure 2).

### Application: Effects of Social Interactions on Participation in Extracurricular Activities
- `00_Add_Health_data.R` prepares the data set to be used.
- `A_Add_Health_nofixed.R` replicates the estimations without fixed effects (Models 1-3).
- `B_Add_Health_fixed.R`replicates the estimations without fixed effects (Models 4-6).
- `C_Add_Health_hete.R` replicates the estimations with heterogeneity in peer effects without addressing network endogeneity (Model 7).
- `D1_Add_Health_Net.BRE.R` replicates the Bayesian estimation of the dyadic linking model (Online Appendix B.3.3).
- `D2_Add_Health_Net.FE.R` replicates the logit estimation of the dyadic linking model.
- `D3_Add_Health_Endo-R` replicates the estimations with heterogeneity in peer effects by addressing network endogeneity (Models 8-9).
- `E_counterfactual.R` replicates the counterfactual analysis.
