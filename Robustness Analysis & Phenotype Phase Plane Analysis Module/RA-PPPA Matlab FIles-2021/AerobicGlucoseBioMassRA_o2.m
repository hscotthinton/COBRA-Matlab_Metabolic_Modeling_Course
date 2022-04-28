% AerobicGlucoseBioMassRA_o2.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set oxygen uptake rate
model = changeRxnBounds(model,'EX_o2(e)',-60,'l'); % set EX_o2(e) = -17

% Set the upper bound for glucose uptake
model = changeRxnBounds(model,'EX_glc(e)',-5,'l'); % set EX_glc(e) = -18.5

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Using robustnessAnalysis, plot the objective function as a function of 
% the oxygen uptake rate
[controlFlux, objFlux] = robustnessAnalysis(model,'EX_o2(e)',100);
