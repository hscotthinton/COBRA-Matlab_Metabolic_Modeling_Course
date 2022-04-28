% AerobicOxygenBioMassRA_Map.m
clear;
% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set uptake rates for robustness analysis of oxygen
model = changeRxnBounds(model,'EX_glc(e)',-10,'b'); % Fix glucose level
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

FBAsolution = optimizeCbModel(model,'max',0,0);

% Import E.coli core map and adjust parameters
map=readCbMap('ecoli_Textbook_ExportMap');
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;
drawFlux(map, model, FBAsolution.x, options);

% Print Fluxes
printFluxVector(model, FBAsolution.x, true)

