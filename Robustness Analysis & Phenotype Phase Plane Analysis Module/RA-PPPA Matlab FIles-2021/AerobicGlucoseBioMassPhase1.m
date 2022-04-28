% AerobicGlucoseBioMassPhase1.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set oxygen and glucose uptake rates
model = changeRxnBounds(model,'EX_glc(e)',-1,'b');
model = changeRxnBounds(model,'EX_o2(e)',-10,'b');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Perform FBA with Biomass_Ecoli_core_N(w/GAM)_Nmet2 as the objective, 
FBAsolution = optimizeCbModel(model,'max',0,1)

% Import E.coli core map and adjust parameters
map=readCbMap('ecoli_Textbook_ExportMap');
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;

% Draw the flux values on the map "target.svg" which can be opened in FireFox
drawFlux(map, model, FBAsolution.x, options);

% Print flux values
printFluxVector(model, FBAsolution.x, true)

