% ShadowPricesAerobicGrowthRate_glc.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set the upper bound for oxygen uptake
model = changeRxnBounds(model,'EX_o2(e)',-20,'b');

% Set the upper bound for glucose uptake
model = changeRxnBounds(model,'EX_glc(e)',-11,'b');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Optimize objective function
FBAsolution = optimizeCbModel(model,'max');
% FBAsolution = optimizeCbModel_old(model,'max'); % Current version does not provide all shadow prices

% Print flux values
printFluxVector(model, FBAsolution.x, true)

% Print shadow prices
'Shadow prices'
printShadowPriceVector(model, - FBAsolution.y, true)

