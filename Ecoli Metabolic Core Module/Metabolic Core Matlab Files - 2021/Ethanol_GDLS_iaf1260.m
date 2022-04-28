% Ethanol_GDLS_iaf1260.m

clear;

% Input the E.coli core model

model=readCbModel('iAF1260.mat');

% Set carbon source and oxygen uptake rates

model = changeRxnBounds(model,'EX_glc(e)',-5,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

% Choose selected reactions
[transRxns,nonTransRxns] = findTransRxns(model,true);
[tmp,ATPMnumber] = ismember('ATPM',nonTransRxns);
[tmp,BioMassnumber] = ismember('Ec_biomass_iAF1260_core_59p81M',nonTransRxns);
nonTransRxnsLength = length(nonTransRxns);
selectedRxns = {nonTransRxns{[1:ATPMnumber-1, ATPMnumber+1:nonTransRxnsLength]}};

% GDLS analysis for Ethanol secretion

[gdlsSolution, bilevelMILPproblem, gdlsSolutionStructs] = GDLS(model, 'EX_etoh(e)', 'minGrowth', 0.03, 'selectedRxns', selectedRxns, 'maxKO', 5, 'nbhdsz', 1);

