% EthanolProduction_multiProductionEnvelope_GDLS_genes.m
clear; clc;

model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% 3 Knockout genes
geneList = {'b1761','b1676','b2297'};
% [model,hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(model,geneList);

[results, ListResults] = findRxnsFromGenes(model, geneList,0,1);

biomassRxn = {'Biomass_Ecoli_core_N(w/GAM)-Nmet2'};
deletions = ListResults(1:3); 

%Show only growth coupled metabolites
figure(1)
[biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,20,false);
% [biomassValues,targetValues] = ProductionEnvelope(model,deletions,biomassRxn,false,20,false);

%Show all secreted metabolites
figure(2)
[biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,20,true);