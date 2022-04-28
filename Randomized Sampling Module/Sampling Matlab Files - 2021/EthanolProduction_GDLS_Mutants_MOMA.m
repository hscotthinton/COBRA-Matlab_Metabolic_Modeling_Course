% EthanolProduction_GDLS_Mutants_MOMA.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set carbon source and oxygen uptake rates for wild type model
model = changeRxnBounds(model,'EX_glc(e)',-5,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');
FBAsolution = optimizeCbModel(model,'max',0,0);
modelWT = model;

% Knockout reactions for mutant model
model = changeRxnBounds(model,'NADH16',0,'b');
model = changeRxnBounds(model,'PTAr',0,'b');
model = changeRxnBounds(model,'TKT2',0,'b');
Mutantsolution = optimizeCbModel(model,'max',0,0);
modelMutant = model;

[solutionDel,solutionWT,totalFluxDiff,solStatus] = MOMA(modelWT,modelMutant,'max',false)

printFluxVector(model, [FBAsolution.x,Mutantsolution.x solutionDel.x], true)
%printFluxVector(model, [FBAsolution.x,solutionWT.x solutionDel.x], true)

% Import E.coli core map and adjust parameters
map=readCbMap('ecoli_Textbook_ExportMap');
options.lb = -10;
options.ub = 10;
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;

% Draw the MOMA flux values on the map "target.svg" which can be opened in FireFox
drawFlux(map, model, solutionDel.x, options);
