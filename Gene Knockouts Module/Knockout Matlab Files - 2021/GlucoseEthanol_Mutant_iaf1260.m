% GlucoseEthanol_Mutant_iaf1260.m
clear; 

% Set constraints
model = readCbModel('iAF1260.mat');
model = changeRxnBounds(model, {'EX_o2(e)', 'EX_glc(e)'}, [0 -20], 'l');
model = changeObjective(model,'Ec_biomass_iAF1260_core_59p81M');

% 1 Knockout reactions
% model = changeRxnBounds(model,'G6PDH2r',0,'b');

% 2 Knockout reactions
% model = changeRxnBounds(model,'ALCD2x',0,'b');
% model = changeRxnBounds(model,'PTAr',0,'b');

% 3 Knockout reactions
model = changeRxnBounds(model,'GLUDy',0,'b');
model = changeRxnBounds(model,'PFL',0,'b');
model = changeRxnBounds(model,'TPI',0,'b');

% Perform FBA  
FBAsolution = optimizeCbModel(model,'max',0,0)
printFluxVector(model,FBAsolution.x, true,true)

% Import E.coli core map and adjust parameters
map=readCbMap('ecoli_iaf1260_ExportMap_Revised');
options.lb = -10;
options.ub = 10;
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;
drawFlux(map, model, FBAsolution.x, options);
