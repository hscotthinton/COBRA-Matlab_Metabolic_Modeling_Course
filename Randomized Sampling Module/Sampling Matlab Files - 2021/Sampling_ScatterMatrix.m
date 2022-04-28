% Sampling_ScatterMatrix.m
clear;

% Input the E.coli core model and set constraints
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
% model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Sample model
sampleStruct = gpSampler(model,5000,[],120);

% Plot scatter matrix
rxnList = {'PGI', 'PFK', 'FBP', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'};
sampleScatterMatrix(rxnList,model,sampleStruct.points,250,24,true);
