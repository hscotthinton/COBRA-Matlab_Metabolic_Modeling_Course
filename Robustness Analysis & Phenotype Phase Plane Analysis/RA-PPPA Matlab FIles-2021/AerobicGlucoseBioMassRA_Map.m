% AerobicGlucoseBioMassRA_Map.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set uptake rates for robustenss analysis of glucose
model = changeRxnBounds(model,'EX_glc(e)',-7,'b'); % Values include -7 and -10
model = changeRxnBounds(model,'EX_o2(e)',-18.5,'l'); % Fix oxygen level

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)_Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

FBAsolution = optimizeCbModel(model,'max',0,0);

% Import E.coli core map and adjust parameters
map=readCbMap('ecoli_Textbook_ExportMap');
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;
drawFlux(map, model, FBAsolution.x, options);

% Print Fluxes
printFluxVector(model, FBAsolution.x, true)

