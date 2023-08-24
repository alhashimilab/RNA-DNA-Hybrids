# RNA-DNA-Hybrids
G-T/U mismatch dynamics in RNA-DNA hybrids
#########################################################################
# RNA-DNA Hybrids
#  Processed Data and Scripts to analyze and plot data for RNA-DNA hybrid manuscript
#  titled:
#  NMR measurements of transient low-populated tautomeric and anionic Watson-Crick-like
#  G•T/U in RNA:DNA hybrids: Implications for the fidelity of transcription and
#  CRISPR/Cas9 gene editing
#  Or Szekely
#  2023-08-23
#########################################################################

This repository contains data folders and scripts for analysis and plotting.

1. 2023-RD-Figures (NMR R1rho data - Figure 2 and Figure S2)
1.1. Raw R1rho data was processed using Disprun (https://github.com/alhashimilab/Disprun.git) based on NMRPipe (Version 10.6 Revision 2020.007.13.34) to extract peak intensities and fit the relaxation curves using mono-exponential decays. The resulting R1rho values are saved in 2023-R1rho-tables.
1.2. R1rho values are then fit to Bloch-McConnell equations using BMNS (https://github.com/alhashimilab/BMNS.git). The data from two probes are either fit individually, or globally, to a 2-state / 3-state exchange model. Resulting fits are saved as csv files in folders under 2023-RD-Figures/2023-RD-RawData.
1.2. Use the script 2023-RD-Figures/plot_hybrid_RD_fits_all.py (python 2.7) to plot the R1rho profiles (R2eff vs. omega_obs) for each one of the cases.
1.3. Use BMNS.py -compare model_1.csv model_2.csv (from BMNS (https://github.com/alhashimilab/BMNS.git)) to compare two models and get AIC and BIC statistical weights for Figure S2.

2. Exchange parameters (Hybrids - Figure 2E, DNA and RNA - Figure S5)
2.1. Exchange parameters for all constructs are saved under 2023-Figures-Main-and-SI/2023-NMR-RD-params.xlsx
2.2. The Matlab script 2023-Figures-Main-and-SI/Figures_all.m will read the exchange parameters and plot bar graphs for Figure 2E and Figure S5 as separate panels.

3. Kinetic rate constants pH extrapolation (Figure S3)
3.1. Use 2023-Figures-Main-and-SI/Figures_all.m to read kinetic rate constants from previously measured pH-dependent DNA data
3.2. Use 2023-Figures-Main-and-SI/Figures_all.m to plot ln(k) vs. pH (Figure S3) and use the Matlab function fit to get a linear fit and extract a and b values (ln(k) = a*pH + b). 

4. Flux calculations at different pHs (Figure 4; Figure S6; Figure S7)
4.1. Within 2023-Figures-Main-and-SI/Figures_all.m call a function (pH_interpolation.m) to use the slope a (from section 3.2) to extrapolate the kinetic rate constants to other pHs.
4.2. Write csv files with the extrapolated kinetic rate constants, save the files under 2023-Figures-Main-and-SI/params_for_python_flux.
4.3. Use a python2.7 script (2023-Figures-Main-and-SI/Flux-Simulations-Pol-Epsilon/flux-simulations-pol-epsilon-k2-vary-all.py) to read the csv file, solve the kinetic differential equation system for the MIS scheme with varying k2s, and calculate the fractional flux through ES2 (fA) for each case.
4.4. Use 2023-Figures-Main-and-SI/Figures_all.m to read the csv files and plot fA vs. k2 for the extrapolated pHs - pH 7.4 (Figure 4) and pH 6.9 (Figure S7A), and the original pHs (Figure S7B).
4.5 Use a python2.7 script (2023-Figures-Main-and-SI/Flux-Simulations-Pol-Epsilon/flux-simulations-pol-epsilon-k2-vary-all.py) to read the csv files for DNA and RNA at pH 8.4, solve the kinetic differential equation system for the MIS scheme with varying k2s and varying k-1, and calculate the fractional flux through ES2 (fA) for each case.
4.6. Use 2023-Figures-Main-and-SI/Figures_all.m to read the csv files and plot an fA heat map for different k2 and k-1 (Figure S6B).

5. NMR Spectra (Figure S1)
5.1. Process 1D and 2D NMR spectra using NMRPipe (Version 10.6 Revision 2020.007.13.34) and convert spectrum to Matlab format (.mat).
5.2. Spectra .mat files are saved in 2023-NMR-Spectra/2023-1D-2D-imino-RawData
5.3. Use 2023-Figures-Main-and-SI/Figures_all.m to load the .mat files and plot 1D and 2D spectra (Figure S1)


6. NTP NMR titrations (Figure 3)
6.1.  Process 2D NMR titration spectra using NMRPipe (Version 10.6 Revision 2020.007.13.34) and convert spectrum to SPARKY format (.ucsf).
6.2. Open ucsf files using SPARKY software. SPARKY projects are saved in /Users/orsula1/OrsDocs/OrsDocs_since_Jan2023/2023-RNA-DNA-Hybrid/2023-Manuscript-Figures/2023-For-Data-Deposition/2023-NMR-Spectra/2023-NMR-titrations/NMR-titrations-Processed-Spectra. Use SPARKY to create figure of overplayed spectra (Figure S4).
6.3. Use SPARKY to pick peaks and save them as SPARKY list files.
6.4. Create .mat files for the chemical shifts of the peaks (saved in 2023-NMR-Spectra/2023-NMR-titrations/cs_dTTP.mat and 2023-NMR-Spectra/2023-NMR-titrations/cs_rUTP.mat)
6.5. Use 2023-Figures-Main-and-SI/Figures_all.m to read the chemical shifts and read the measured pH values (saved in 2023-NMR-Spectra/2023-NMR-titrations/dTTP-rUTP-actual-pH.xlsx). Plot chemical shift change vs. pH, and use Matlab function nlinfit to fit the data to the chemical shift change equation (Eq. 5) and extract pKa.
