% AerobicSuccinateBioMass.m

clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set maximum uptake rate for carbon sources:

model = changeRxnBounds(model,'EX_glc(e)',0,'l');
model = changeRxnBounds(model,'EX_succ(e)',-20,'l');

% Set oxygen uptake rate
model = changeRxnBounds(model,'EX_o2(e)',-1000,'l');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)_Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Perform FBA with Biomass_Ecoli_core_N(w/GAM)_Nmet2 as the objective, 
FBAsolution = optimizeCbModel(model,'max',0,0)

% Print fluxes
printFluxVector(model, FBAsolution.x, true)

% Import E.coli core map and adjust parameters
map=readCbMap('ecoli_Textbook_ExportMap');
options.lb = -10;
options.ub = 10;
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;

% Draw the flux values on the map "target.svg" which can be opened in FireFox
drawFlux(map, model, FBAsolution.x, options);
