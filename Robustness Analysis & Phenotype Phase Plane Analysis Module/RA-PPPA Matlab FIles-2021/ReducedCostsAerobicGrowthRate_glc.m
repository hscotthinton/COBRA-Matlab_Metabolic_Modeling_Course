% ReducedCostAerobicGrowthRate_glc.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set the lower bounds for oxygen adn glucose uptake
model = changeRxnBounds(model,'EX_o2(e)',-20,'b');
model = changeRxnBounds(model,'EX_glc(e)',-10,'b');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Optimize objective function
FBAsolution = optimizeCbModel(model,'max');

% Print flux values

'Flux values and Reduced Costs'
printFluxVector(model, [FBAsolution.x,FBAsolution.w], true)


