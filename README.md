# sEMG Muscle Activation Pipeline
MATLAB scripts for processing surface electromyography (sEMG) data from cylindrical grasp trials (10 participants, 5 trials each) for muscle activation analysis supporting the paper at https://doi.org/10.1016/j.jbiomech.2025.112060.​

# Pipeline Overview
The workflow compares "Old" (windowed envelope) vs "New" (Hilbert-Huang Transform + ODE) pipelines for EDC, FDS, FDP muscles, then runs simplified linear mixed-effects (LME) models per muscle: log(activation) ~ pipeline + trial + (1|participant).​

# File Descriptions
Script	Purpose	Key Outputs
sEMG_2_Activations.m	Main processing: HHT filtering, ODE activation, Excel export	activations_export_marginal.xlsx, plots (activations, ratios, HHT spectra), Appendix CSV
activation.m	ODE model from Blana et al. (2007): da/dt = (u/Tact + (1-u)/Tdeact)(u - a)	Called by main script
LME_muscle_activations_long_file_creation.m	Converts Excel to LME-ready long CSV (handles skips like P5-T5)	activations_long_raw_marginal.csv
run_lme_activation_simple.m	Runs per-muscle LME (REML), winsorizes log(activation), bins time	LME CSVs (fixed effects, ratios, residuals plots), p-value adjustments
Usage Instructions
Place sEMG Excel files (sEMG_P1.xlsx to sEMG_P10.xlsx) in the data path (updates to your OneDrive/University of Warwick/Hand Trials/Results/Cylindrical Grasp/sEMG Data).​

Run sEMG_2_Activations.m first to generate activations_export_marginal.xlsx and plots.

Run LME_muscle_activations_long_file_creation.m to create activations_long_raw_marginal.csv.

Run run_lme_activation_simple.m for LME results in LME_Results/ (creates folders automatically).​

Adjust Zero=1 (zero initial condition), HHT_View=1 (enable HHT plots), or NBINS_PER_TRIAL=20 as needed.​

# Dependencies & Notes
MATLAB Signal Processing Toolbox (hht, emd, butter, filtfilt), Statistics Toolbox (fitlme), Curve Fitting Toolbox (spline, mkpp, ppval).

Results include residual diagnostics, Holm-Bonferroni p-values, and geometric mean ratios (New/Old).​
