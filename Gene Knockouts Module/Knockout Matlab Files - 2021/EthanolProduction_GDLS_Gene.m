% EthanolProduction_GDLS_Gene.m
clear; 

% Set operating conditions
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',0,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Select reactions
[transRxns,nonTransRxns] = findTransRxns(model,true);
[tmp,ATPMnumber] = ismember('ATPM',nonTransRxns); % Identify ATPM reaction number
[tmp,BioMassnumber] = ismember('Biomass_Ecoli_core_N(w/GAM)-Nmet2',nonTransRxns); % Identify biomass reaction number
nonTransRxnsLength = length(nonTransRxns); % Find number of non-transport reactions
selectedRxns = {nonTransRxns{[1:ATPMnumber-1, ATPMnumber+1:BioMassnumber-1, BioMassnumber+1:nonTransRxnsLength]}}; 

% [geneList]=findGenesFromRxns(model,selectedRxns)
selectedRxns = model.genes;

% GDLS analysis for Ethanol secretion
[gdlsSolution, bilevelMILPproblem, gdlsSolutionStructs] = GDLS(model, 'EX_etoh(e)', 'minGrowth', 0.05, ...
    'koType','genes','selectedRxns', selectedRxns, 'maxKO', 3, 'nbhdsz', 2);

gdlsSolution.KOs
[results ListResults] = findRxnsFromGenes(model, gdlsSolution.KOs)


% options.selectedRxns = selectedRxns;
% options.minGrowth = 0.05;
% koType = 'rxns';
% % selectedRxns = model.genes;
% options.maxKO = 1;
% options.nbhdsz = 2;
% [gdlsSolution, bilevelMILPproblem, gdlsSolutionStructs] = GDLS(model,'EX_etoh(e)');
