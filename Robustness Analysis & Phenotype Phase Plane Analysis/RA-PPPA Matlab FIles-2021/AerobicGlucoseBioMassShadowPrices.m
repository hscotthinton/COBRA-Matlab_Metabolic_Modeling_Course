% AerobicGlucoseBioMassShadowPrices.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

model = changeRxnBounds(model,'EX_glc(e)',-10,'b');
model = changeRxnBounds(model,'EX_o2(e)',-20,'b'); 

% Set optimization objective
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Calculate flux values
FBAsolution = optimizeCbModel(model,'max') % You can't remove loops and calculate shadow prices at the same time!
% FBAsolution = optimizeCbModel_old(model,'max'); % Current version does not provide all shadow prices

% Print flux values and reduced costs
disp('  ')
disp('Flux values and Reduced Costs')
printFluxVector(model, [FBAsolution.x,FBAsolution.w], true)

% Print Shadow Prices
disp('  ')
disp('Shadow Prices')
% printShadowPriceVector(model, FBAsolution.y, true) % Use for Tomlab_CPLEX solvers
printShadowPriceVector(model, - FBAsolution.y, true) % Use for gurobi and gplk solver

% Import E.coli core map and adjust parameters
map=readCbMap('ecoli_Textbook_ExportMap');
options.lb = -10;
options.ub = 10;
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;
drawFlux(map, model, FBAsolution.x, options);