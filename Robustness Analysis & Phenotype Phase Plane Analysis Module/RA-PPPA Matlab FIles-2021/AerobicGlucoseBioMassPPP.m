% AerobicGlucoseBioMassPPP.m

clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set the lower bounds for oxygen and glucose uptake
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
model = changeRxnBounds(model,'EX_glc(e)',-20,'l');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Using phenotype plane analysis, plot the objective function as a 
% function of the glucose and oxygen uptake rate
[growthRates,shadowPrices1,shadowPrices2]=phenotypePhasePlane(model,'EX_glc(e)','EX_o2(e)');
