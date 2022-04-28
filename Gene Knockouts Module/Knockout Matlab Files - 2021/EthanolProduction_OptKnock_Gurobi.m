% EthanolProduction_OptKnock_gurobi.m
clear;

% Load the E.coli core model
model = readCbModel('ecoli_core_model.mat');
solverOK = changeCobraSolver('gurobi','all');

% Set carbon source and oxygen uptake rates
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

% Set optimization objective
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Select reactions to be explored by the optKnock algorithm
[transRxns,nonTransRxns] = findTransRxns(model,true);
[tmp,ATPMnumber] = ismember('ATPM',nonTransRxns); % Identify ATPM reaction number
[tmp,BioMassnumber] = ismember('Biomass_Ecoli_core_N(w/GAM)-Nmet2',nonTransRxns); % Identify biomass reaction number
nonTransRxnsLength = length(nonTransRxns); % Find number of non-transport reactions
selectedRxns = {nonTransRxns{[1:ATPMnumber-1, ATPMnumber+1:BioMassnumber-1, BioMassnumber+1:nonTransRxnsLength]}}; 

% Optknock analysis for Ethanol secretion
options.targetRxn = 'EX_etoh(e)';
options.vMax = 1000;
options.numDel = 3;
options.numDelSense = 'L';
constrOpt.rxnList = {'Biomass_Ecoli_core_N(w/GAM)-Nmet2','ATPM'};
constrOpt.values = [0.14, 8.39]; % Changing minimum growth-rate (biomass) from 0.05 to 0.14
constrOpt.sense = 'GE';

optKnockSol = OptKnock(model, selectedRxns, options, constrOpt);
deletions = optKnockSol.rxnList'

% Print out growth rate and minimum & maximum secretion rate
[growthRate,minProd,maxProd] = testOptKnockSol(model,'EX_etoh(e)',optKnockSol.rxnList)

% Print production envelope
lineColor = 'b';
targetRxn = 'EX_etoh(e)';
biomassRxn = 'Biomass_Ecoli_core_N(w/GAM)-Nmet2';
geneDelFlag = false;
nPts = 50;

[biomassValues,targetValues] = productionEnvelope(model,deletions,lineColor,targetRxn,biomassRxn,geneDelFlag,nPts);
xlabel('Biomass (mmol/g DW-hr)')
ylabel('EX-etoh(e)(mmol/g DW-hr)')

