% Sampling_sampleCbModel_Example_gpSampler.m
clear;clc;

% Input the E.coli core model and set constraints
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
biomassRxnAbbr = 'Biomass_Ecoli_core_N(w/GAM)-Nmet2';
ibm = find(ismember(model.rxns, biomassRxnAbbr));  % column index of the biomass reaction
model.lb(ibm)=0.05;
model.c(:)=0; % Remove biomass as objective function

model_aerobic = changeRxnBounds(model,'EX_o2(e)',-20,'l');
model_anaerobic = changeRxnBounds(model_aerobic,'EX_o2(e)',0,'l');

options.toRound = 1;
options.nFiles = 50;
options.nPointsReturned = 5000;
[aerobic_modelSampling,aerobic_samples] = sampleCbModel(model_aerobic,[],'ACHR',options);
[anaerobic_modelSampling,anaerobic_samples] = sampleCbModel(model_anaerobic,[],[],options);

% Visualize sampling results for a set of reactions. 
rxnList = {'PGI', 'PFK', 'FBP', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'};
plotSampleHist(rxnList, {aerobic_samples, anaerobic_samples}, {model_aerobic, model_anaerobic},[],[5,2]);