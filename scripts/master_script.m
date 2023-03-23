
%% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% TB preventive therapy in South Africa
% All code writen by Juan F Vesga (juan.vesga-gaviria@lshtm.ac.uk)
%% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

%% Necesary first steps
clear
close all

% Add folder and subfolders to path
addpath(genpath(pwd));

%Call Global parameters necesary to build models and retrieve data data
get_globals; % Explore "get_globals.m" for more details

%% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% 1)  Plot model calibrations for South Africa against WHO data

% Load runs produced by sampling the posterior after MCMC calibration 
file=sprintf('%s','output','_',location,'_','mcmc','.mat'); 
load(file);
runs=object;

% Run this line to plot model vs data and save graph in "results"
plot_targets(runs,data,location); % calls a plot function in "functions" folder


%% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% 2)  Plot pre-run analysis of Partial Rank Correlation of regimen
% attributes

% Select scenarios 
scenario = "base"; % options: "base" and "rif_free" (fo rifampin free regimens)
PRCC_plots;









