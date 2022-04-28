% GlucoseEthanolOptknock_iaf1260_Reduced_Reactions.m
clear; 

% Set constraints
model = readCbModel('iAF1260.mat');
model = changeRxnBounds(model, {'EX_o2(e)', 'EX_glc(e)'}, [0 -10], 'l');
model = changeObjective(model,'Ec_biomass_iAF1260_core_59p81M');

[transRxns,nonTransRxns] = findTransRxns(model,true); % Remove transport reactions
includedSubSystems = {'Transport, Inner Membrane','Glycerophospholipid Metabolism','Transport, Outer Membrane Porin',... 
    'Cell Envelope Biosynthesis','Nucleotide Salvage Pathway','Murein Recycling','Membrane Lipid Metabolism', ...
    'Glycerophospholipid Metabolism','Inorganic Ion Transport and Metabolism','Lipopolysaccharide Biosynthesis / Recycling',...
    'tRNA Charging','Unassigned','Membrane Lipid Metabolism','Murein Biosynthesis'}; 
% unwantedReactions = model.rxns(ismember(model.subSystems,includedSubSystems));
[unwantedReactions,rxnPos] = findRxnsFromSubSystem(model,includedSubSystems)
[tf,rids] = ismember([unwantedReactions;{'ATPM'};{'Ec_biomass_iAF1260_core_59p81M'}], nonTransRxns);
idfinal=rids(tf);
nonTransRxns(idfinal) = [];
selectedRxns = nonTransRxns; 

% Optknock analysis for Ethanol secretion
disp('Executing optKnock');
options.targetRxn = 'EX_etoh(e)';
options.vMax = 1000;
options.numDel = 3;
options.numDelSense = 'L';
constrOpt.rxnList = {'Ec_biomass_iAF1260_core_59p81M','ATPM'};
constrOpt.values = [0.05, 8.39];
constrOpt.sense = 'GE';

optKnockSol = OptKnock(model, selectedRxns, options, constrOpt);
deletions = optKnockSol.rxnList'

% Print out growth rate and minimum & maximum secretion rate
[growthRate,minProd,maxProd] = testOptKnockSol(model,'EX_etoh(e)',optKnockSol.rxnList)

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

