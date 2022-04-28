% AnaerobicEthanolRA.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set uptake rates
model = changeRxnBounds(model,'EX_glc(e)',-10,'b'); 
model = changeRxnBounds(model,'EX_o2(e)',-0,'b'); 
% model = changeRxnBounds(model,'EX_etoh(e)',0,'l'); 

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Using robustnessAnalysis, plot the objective function as a function 
% of the ethanol secretion rate

[controlFlux, objFlux] = robustnessAnalysis(model,'EX_etoh(e)',100);
