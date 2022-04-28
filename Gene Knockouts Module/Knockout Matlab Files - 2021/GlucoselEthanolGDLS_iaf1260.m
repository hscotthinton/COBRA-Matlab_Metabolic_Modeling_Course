% GlucoseEthanolGDLS_iaf1260.m
clear; clc;

% Set constraints
model = readCbModel('iAF1260.mat');

model = changeRxnBounds(model, {'EX_o2(e)', 'EX_glc(e)'}, [0 -10], 'l');
model = changeObjective(model,'Ec_biomass_iAF1260_core_59p81M');

% Remove exchange & transport reactions & ATPM & Biomass reaction
[transRxns,nonTransRxns] = findTransRxns(model,true); % Remove transport reactions
[tmp,ATPMnumber] = ismember('ATPM',nonTransRxns); % Identify ATPM reaction number
nonTransRxns(ATPMnumber) = []; % Remove ATPM from potential knockout reactions
selectedRxns = nonTransRxns;

% GDLS analysis for Ethanol secretion
[gdlsSolution, bilevelMILPproblem, gdlsSolutionStructs] = GDLS(model, 'EX_etoh(e)', 'minGrowth', 0.05, ...
    'selectedRxns', selectedRxns, 'maxKO', 3, 'nbhdsz', 2);

[growthRate,minProd,maxProd] = testOptKnockSol(model,'EX_etoh(e)',gdlsSolution.KOs)
deletions = gdlsSolution.KOs'

% Print production envelope
lineColor = 'b';
targetRxn = 'EX_etoh(e)';
biomassRxn = 'Ec_biomass_iAF1260_core_59p81M';
geneDelFlag = false;
nPts = 50;
figure(1)
[biomassValues,targetValues] = productionEnvelope(model,deletions,lineColor,targetRxn,biomassRxn,geneDelFlag,nPts);
xlabel('Biomass (mmol/g DW-hr)')
ylabel('Ethanol(mmol/g DW-hr)')

figure(2)
[biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,50,false);

