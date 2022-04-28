% EthanolProduction_GDLS_Gene_iaf1260.m
clear; 

% Set operating conditions
model = readCbModel('iAF1260.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',0,'l');
model = changeObjective(model,'Ec_biomass_iAF1260_core_59p81M');

selectedRxns = model.genes;

% GDLS analysis for Ethanol secretion
[gdlsSolution, bilevelMILPproblem, gdlsSolutionStructs] = GDLS(model, 'EX_etoh(e)', 'minGrowth', 0.05, ...
    'koType','genes','selectedRxns', selectedRxns, 'maxKO', 3, 'nbhdsz', 2);

gdlsSolution.KOs
[results ListResults] = findRxnsFromGenes(model, gdlsSolution.KOs)


