% AerobicGlucoseBioMassPhaseOptimal.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set oxygen and glucose uptake rates
model = changeRxnBounds(model,'EX_glc(e)',-3,'b');
model = changeRxnBounds(model,'EX_o2(e)',-10,'l');

% Preventing fermentation
model = changeRxnBounds(model,'EX_ac(e)',-0,'b');
model = changeRxnBounds(model,'EX_for(e)',-0,'b');
model = changeRxnBounds(model,'EX_etoh(e)',-0,'b');
model = changeRxnBounds(model,'EX_lac-D(e)',-0,'b');
model = changeRxnBounds(model,'EX_pyr(e)',-0,'b');
model = changeRxnBounds(model,'EX_acald(e)',-0,'b');

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

