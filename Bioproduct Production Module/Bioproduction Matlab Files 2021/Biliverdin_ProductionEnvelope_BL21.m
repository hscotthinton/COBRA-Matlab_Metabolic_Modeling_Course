% Biliverdin_ProductionEnvelope_BL21.m
clear;

% Input the E.coli core model
model=readCbModel('iECD_1391.xml');

% Add heme oxygenase enzyme
model=addReaction(model,'HEMEOX','pheme[c] + 3 nadph[c] + 5 h[c] + 3 o2[c] -> biliverdin[c] + fe2[c] + co[c] + 3 nadp[c] + 3 h2o[c]');

% Add biliverdin secretion reactions
model = addReaction(model,'BILIVERDINtex','biliverdin[p] <=> biliverdin[e]');
model = addReaction(model,'BILIVERDINtpp','biliverdin[c] <=> biliverdin[p]');
model = addReaction(model,'EX_biliverdin(e)','biliverdin[e] <=>');

metID=findMetIDs(model,{'biliverdin[c]', 'biliverdin[p]', 'biliverdin[e]'});
model.metCharge(metID(1))= -2;
model.metCharge(metID(2))= -2;
model.metCharge(metID(3))= -2;

% Add carbon monoxide secretion reactions
model = addReaction(model,'COtpp','co[c] <=> co[p]');
model = addReaction(model,'COtex','co[p] <=> co[e]');
model = addReaction(model,'EX_co(e)','co[e] <=>');
rxnID=findRxnIDs(model,'COtpp');
model.confidenceScores(rxnID)={0};
rxnID=findRxnIDs(model,'COtex');
model.confidenceScores(rxnID)={0};

% Define confidence levels
rxnID=findRxnIDs(model,'HEMEOX');
model.confidenceScores(rxnID)={0};
rxnID=findRxnIDs(model,'BILIVERDINtex');
model.confidenceScores(rxnID)={0};
rxnID=findRxnIDs(model,'BILIVERDINtpp');
model.confidenceScores(rxnID)={0};

% Set key variables
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
model = changeObjective(model,'Ec_biomass_iJO1366_core_53p95M');
% model = changeRxnBounds(model,'EX_biliverdin(e)',0.5,'l');

deletions = {}; 
lineColor = 'b';
targetRxn = 'EX_biliverdin(e)';
biomassRxn = 'Ec_biomass_iJO1366_core_53p95M';
geneDelFlag = false;
nPts = 50;

[biomassValues,targetValues] = productionEnvelope(model,deletions,lineColor,targetRxn,biomassRxn,geneDelFlag,nPts);
title('Biliverdin Production Envelope'); xlabel('Biomass (mmol/g DW-hr)'); ylabel('EX-biliverdin(e)(mmol/g DW-hr)')
axis tight
