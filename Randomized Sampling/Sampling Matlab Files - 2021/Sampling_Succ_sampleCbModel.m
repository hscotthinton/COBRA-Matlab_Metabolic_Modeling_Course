% Sampling_Succ_sampleCbModel.m
clear; 

model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',0,'l');
model = changeRxnBounds(model,'EX_o2(e)',-40,'l');
model = changeRxnBounds(model,'EX_succ(e)',-20,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% model = changeRxnBounds(model,'Biomass_Ecoli_core_N(w/GAM)_Nmet2',0.84,'b');
FBAsolution = optimizeCbModel(model,'max',0,0);
printFluxVector(model, FBAsolution.x, true);

%1) Sample a model called 'superModel' using default settings and save the
%   results in files with the common beginning 'superModelSamples'

% [modelSampling,samples] = sampleCbModel(model,'superModelSamples');
% [modelSampling,samples] = sampleCbModel(model);

% %2) Sample a model called 'hyperModel' using default settings except with a total of 50 sample files
% %   saved and with 5000 sample points returned.
% 
% options.useFastFVA = 0;

biomassRxnAbbr = 'Biomass_Ecoli_core_N(w/GAM)-Nmet2';
ibm = find(ismember(model.rxns, biomassRxnAbbr));  % column index of the biomass reaction
model.lb(ibm)=0.05;
model.c(:)=0;

options.toRound = 1;
options.nFiles = 50;
options.nPointsReturned = 10000;
[modelSampling,samples] = sampleCbModel(model,[],[],options);

% Visualize sampling results for a set of reactions. 
rxnList = {'PGI', 'PFK', 'FBP', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'};
plotSampleHist(rxnList, {samples}, {model},[],[2,5]);
