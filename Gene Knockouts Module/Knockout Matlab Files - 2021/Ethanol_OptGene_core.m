% Ethanol_OptGene_core.m
clear; clc;
changeCobraSolver('gurobi','all');

% Set conditions
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',0,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Select reactions
[transRxns,nonTransRxns] = findTransRxns(model,true);
[tmp,ATPMnumber] = ismember('ATPM',nonTransRxns); % Identify ATPM reaction number
[tmp,BioMassnumber] = ismember('Biomass_Ecoli_core_N(w/GAM)-Nmet2',nonTransRxns); % Identify biomass reaction number
nonTransRxnsLength = length(nonTransRxns); % Find number of non-transport reactions
selectedRxnList = {nonTransRxns{[1:ATPMnumber-1, ATPMnumber+1:BioMassnumber-1, BioMassnumber+1:nonTransRxnsLength]}}; 

% Run optGene
[~, ~, ~, optGeneSol] = optGene(model, 'EX_etoh(e)', 'EX_glc(e)', selectedRxnList, 'MaxKOs', 10, 'TimeLimit', 5000);
optGeneSol.rxnList

% Graph production envelope
figure(2)
lineColor = 'b';
targetRxn = 'EX_etoh(e)';
biomassRxn = 'Biomass_Ecoli_core_N(w/GAM)-Nmet2';
geneDelFlag = false;
nPts = 50;
[biomassValues,targetValues] = productionEnvelope(model,optGeneSol.rxnList,lineColor,targetRxn,biomassRxn,geneDelFlag,nPts);
xlabel('Biomass (mmol/g DW-hr)')
ylabel('EX-etoh(e)(mmol/g DW-hr)')


