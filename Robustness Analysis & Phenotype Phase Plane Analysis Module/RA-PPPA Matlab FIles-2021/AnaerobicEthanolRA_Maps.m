% AnaerobicEthanolRA_Maps.m 
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set uptake rates
model = changeRxnBounds(model,'EX_glc(e)',-10,'b'); 
model = changeRxnBounds(model,'EX_o2(e)',-0,'b'); 

% Set ethanol production rate
model = changeRxnBounds(model,'EX_etoh(e)',4.242,'b'); % Low production (growth rate = 0.1757)
% model = changeRxnBounds(model,'EX_etoh(e)',12.53,'b'); % Medium production (growth rate = 0.2002)
% model = changeRxnBounds(model,'EX_etoh(e)',18.38,'b'); % High production (growth rate = 0.1063) 

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Perform FBA with Biomass_Ecoli_core_N(w/GAM)_Nmet2 as the objective, 

FBAsolution = optimizeCbModel(model,'max',0,0)

% Import E.coli core map and adjust parameters
map=readCbMap('ecoli_Textbook_ExportMap');
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;

% Draw the flux values on the map "target.svg" which can be opened in FireFox
drawFlux(map, model, FBAsolution.x, options);

% Print Fluxes
printFluxVector(model, FBAsolution.x, true)

