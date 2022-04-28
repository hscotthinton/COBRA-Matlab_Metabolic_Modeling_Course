% findingOptimalSolutionsSuccVisualize.m
clear; clc;

model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',0,'l');
model = changeRxnBounds(model,'EX_o2(e)',-40,'l');
model = changeRxnBounds(model,'EX_succ(e)',-20,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% List optimal solutions
solverOK = changeCobraSolver('glpk','MILP');

[solutions] = enumerateOptimalSolutions(model); 

v = solutions.fluxes(:,1); % Select which vector wanted to be mapped (1-3)
printFluxVector(model, v, true)

map=readCbMap('ecoli_Textbook_ExportMap');
options.lb = -10;
options.ub = 10;
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;
drawFlux(map, model, v, options);


