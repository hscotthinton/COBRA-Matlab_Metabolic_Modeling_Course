% AerobicGlucoseBioMassRA.m
clear;
% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set uptake rates for robustenss analysis of glucose
model = changeRxnBounds(model,'EX_glc(e)',-18.5,'l'); 
model = changeRxnBounds(model,'EX_o2(e)',-17,'l'); % Fix oxygen level

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)_Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Using robustnessAnalysis, plot the objective function as a function of the glucose uptake rate
[controlFlux, objFlux] = robustnessAnalysis(model,'EX_glc(e)',100);
