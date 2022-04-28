% GlucoseEthanolOptknock_iJO1366_Reduced_Reactions.m
clear; 

% Set constraints
model = readCbModel('iJO1366.mat');
model = changeRxnBounds(model, {'EX_o2_e', 'EX_glc__D_e'}, [0 -10], 'l');
model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');

[transRxns,nonTransRxns] = findTransRxns(model,true); % Remove transport reactions
includedSubSystems = {'Transport, Inner Membrane','Glycerophospholipid Metabolism','Transport, Outer Membrane Porin',... 
    'Cell Envelope Biosynthesis','Nucleotide Salvage Pathway','Murein Recycling','Membrane Lipid Metabolism', ...
    'Glycerophospholipid Metabolism','Inorganic Ion Transport and Metabolism','Lipopolysaccharide Biosynthesis / Recycling',...
    'tRNA Charging','Unassigned','Membrane Lipid Metabolism','Murein Biosynthesis'}; 
% unwantedReactions = model.rxns(ismember(model.subSystems,includedSubSystems));
[unwantedReactions,rxnPos] = findRxnsFromSubSystem(model,includedSubSystems);
[tf,rids] = ismember([unwantedReactions;{'ATPM'};{'Ec_biomass_iAF1260_core_59p81M'}], nonTransRxns);
idfinal=rids(tf);
nonTransRxns(idfinal) = [];
selectedRxns = nonTransRxns; 

% Optknock analysis for Ethanol secretion
disp('Executing optKnock');
options.targetRxn = 'EX_etoh_e';
options.vMax = 1000;
options.numDel = 3;
options.numDelSense = 'L';
constrOpt.rxnList = {'BIOMASS_Ec_iJO1366_core_53p95M','ATPM'};
constrOpt.values = [0.05, 8.39];
constrOpt.sense = 'GE';

optKnockSol = OptKnock(model, selectedRxns, options, constrOpt);
deletions = optKnockSol.rxnList'

% Print out growth rate and minimum & maximum secretion rate
[growthRate,minProd,maxProd] = testOptKnockSol(model,'EX_etoh_e',optKnockSol.rxnList)

% Print production envelope
lineColor = 'b';
targetRxn = 'EX_etoh_e';
biomassRxn = 'BIOMASS_Ec_iJO1366_core_53p95M';
geneDelFlag = false;
nPts = 50;
figure(1)
[biomassValues,targetValues] = productionEnvelope(model,deletions,lineColor,targetRxn,biomassRxn,geneDelFlag,nPts);
xlabel('Biomass (mmol/g DW-hr)')
ylabel('Ethanol(mmol/g DW-hr)')

figure(2)
[biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,50,false);

