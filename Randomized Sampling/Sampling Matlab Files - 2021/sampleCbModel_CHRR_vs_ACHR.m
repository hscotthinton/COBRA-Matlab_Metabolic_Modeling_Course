% Sampling_sampleCbModel_Example_gpSampler.m
clear;clc;

% Input the E.coli core model and set constraints
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
biomassRxnAbbr = 'Biomass_Ecoli_core_N(w/GAM)-Nmet2';
ibm = find(ismember(model.rxns, biomassRxnAbbr));  % column index of the biomass reaction
model.lb(ibm)=0.05;
model.c(:)=0; % Remove biomass as objective function

options.toRound = 1;
options.nFiles = 50;
options.WarmupPoints = 5000
options.nPointsReturned = 5000;
[modelSampling_ACHR,ACHR_samples] = sampleCbModel(model,[],'ACHR',options);
[modelSampling_CHRR,CHRR_samples] = sampleCbModel(model,[],[],options);

% Visualize sampling results for a set of reactions. 
rxnList = {'PGI', 'PFK', 'FBP', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'};
plotSampleHist(rxnList, {ACHR_samples, CHRR_samples}, {modelSampling_ACHR, model},[],[5,2]);