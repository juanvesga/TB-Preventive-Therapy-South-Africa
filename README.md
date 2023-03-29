# TB-Preventive-Therapy-South-Africa
TB transmission model exploring target product profiles of TPT of future regimens, applied to ZAF

# Background
This MATLAB code describes the transmission of two strains (Drug resistant & drug sensitive) of M.tuberculosis at the population level. 
We use a compartmental deterministic model to capture the transmission and natural history of disease, 
as well as the cascade of care (i.e., care seeking, diagonosis, treatment and cure). 

Importantly, this framework was used to explore in detail the importance of different attributes 
of future preventive therapy (TPT) regimens (e.g, efficacy, duration, resistance barrier etc), and the population level 
effect that results from rolling out TPT interventions in high burden settings. 

Here we provide pre-loaded model fits (this repository does not containthe calibration procedures) for the case of South Africa. This 
can be loaded (see steps below) to produce analysis and results. A detailed description of the model and its method can be found in:

Vesga JF, Lienhardt C, Nsengiyumva P, Campbell JR, Oxlade O, den Boon S, Falzon D, Schwartzman K, 
Churchyard G, Arinaminpathy N. Prioritising attributes for tuberculosis preventive treatment regimens: a modelling analysis. 
BMC Med. 2022 May 18;20(1):182. doi: 10.1186/s12916-022-02378-1. PMID: 35581650; PMCID: PMC9115962. 

It also supports results on cost-effectiveness analysis found in: 

Nsengiyumva NP, Campbell JR, Oxlade O, Vesga JF, Lienhardt C, Trajman A, Falzon D, Den Boon S, Arinaminpathy N, Schwartzman K. 
Scaling up target regimens for tuberculosis preventive treatment in Brazil and South Africa: An analysis of costs and cost-effectiveness. 
PLoS Med. 2022 Jun 13;19(6):e1004032. doi: 10.1371/journal.pmed.1004032. PMID: 35696431; PMCID: PMC9239450.

# Instructions

The code is split in scripts and functions. All data and preliminary result structures are here contained. 
First, open the script called "master_script" in the folder "scripts". 
All actions and functions will be executed from this script.
The folder scripts contain different independent scripts that can be run to produce results, as follows:

1) Plot model calibrations to WHO data
2) Plot Partial Rank Correlation analysis
3) Simulate secanrios of TPT roll-out and plot results    


