% AerobicGlucoseBioMassEscher.m

clear; clc;

% Input the E.coli core model
model = readCbModel('e_coli_core.mat');

% Set carbon sources
model = changeRxnBounds(model,'EX_glc__D_e',0,'l');
model = changeRxnBounds(model,'EX_pyr_e',-20,'l');

% Set for anaarobic conditions:
model = changeRxnBounds(model,'EX_o2_e',0,'l');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'BIOMASS_Ecoli_core_w_GAM');

% Perform FBA with Biomass_Ecoli_core_N(w/GAM)_Nmet2 as the objective, 
FBAsolution = optimizeCbModel(model,'max',0,0) 

% Create table of flux values for Escher map
Reactions = model.rxns;
Flux = round(FBAsolution.x,3);
T = table(Reactions,Flux,'RowNames',model.rxns);
% Write to CSV file
writetable(T,'escher_flux.csv');

% Print flux values
printFluxVector(model, FBAsolution.x, true)

