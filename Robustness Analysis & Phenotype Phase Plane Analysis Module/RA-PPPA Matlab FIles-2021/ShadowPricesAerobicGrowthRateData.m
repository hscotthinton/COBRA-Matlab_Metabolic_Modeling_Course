% ShadowPricesAerobicGrowthRateData.m

clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set the lower bounds for oxygen & glucose uptake
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'b'); % Looking at data for EX_o2(e)=20

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Optimize objective function
% FBAsolution = optimizeCbModel(model,'max'); % Cannot remove loops when calcuating shadow prices
FBAsolution = optimizeCbModel_old(model,'max'); % Current version does not provide all shadow prices

% Print flux values
printFluxVector(model, FBAsolution.x, true)

% Print shadow prices
disp(' ')
disp('Shadow prices')
printShadowPriceVector(model, - FBAsolution.y, true) % Include minus sign to have consistent sign for shadow prices
% printShadowPriceVector(model, FBAsolution.y, true) % No minus sign for Tomlab_cplex solver
