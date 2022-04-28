% AerobicOxygenBioMassRA.m
clear;
% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set uptake rates for robustenss analysis of oxygen
model = changeRxnBounds(model,'EX_glc(e)',-10,'b'); % Fix glucose level
model = changeRxnBounds(model,'EX_o2(e)',-25,'l');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)_Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Using robustnessAnalysis, plot the objective function as a function of the glucose uptake rate
[controlFlux, objFlux] = robustnessAnalysis(model,'EX_o2(e)',100)