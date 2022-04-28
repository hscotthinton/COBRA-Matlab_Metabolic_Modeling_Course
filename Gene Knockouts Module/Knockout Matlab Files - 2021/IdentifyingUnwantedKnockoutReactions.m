% IdentifyingUnwantedKnockoutReactions.m
clear; clc;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set carbon source and oxygen uptake rates
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Perform FBA with Biomass_Ecoli_core_N(w/GAM)-Nmet2 as the objective, 
FBAsolution = optimizeCbModel(model,'max')

% Identify non-tranport reactions
[transRxns,nonTransRxns] = findTransRxns(model,true);

% Removing ATPM and biomass function
[tmp,ATPMnumber] = ismember('ATPM',nonTransRxns); % Identify ATPM reaction number
[tmp,BioMassnumber] = ismember('Biomass_Ecoli_core_N(w/GAM)-Nmet2',nonTransRxns); % Identify biomass reaction number
nonTransRxnsLength = length(nonTransRxns); % Find number of non-transport reactions
selectedRxns = {nonTransRxns{[1:ATPMnumber-1, ATPMnumber+1:BioMassnumber-1, BioMassnumber+1:nonTransRxnsLength]}}; 

disp('Transport Reactions');
disp(transRxns);

disp('Nontransport Reactions');
disp(nonTransRxns);

disp('selectedRxns');
disp(selectedRxns');