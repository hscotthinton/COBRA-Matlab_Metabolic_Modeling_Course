% EthanolProduction_multiProductionEnvelope_GDLS_genes_iaf1260.m
clear; clc;

model = readCbModel('iAF1260.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');
model = changeObjective(model,'Ec_biomass_iAF1260_core_59p81M');

% 3 Knockout genes
geneList = {'b1761','b3919','b3951'};
[results, ListResults] = findRxnsFromGenes(model, geneList,0,1);

deletions = ListResults(1:3); 
biomassRxn = {'Ec_biomass_iAF1260_core_59p81M'};

%Show only growth coupled metabolites
figure(1)
[biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,20,false);
% [biomassValues,targetValues] = ProductionEnvelope(model,deletions,biomassRxn,false,20,false);

% %Show all secreted metabolites
% figure(2)
% [biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,20,true);