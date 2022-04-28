% Biliverdin_Secretions_Table_BL21.m
clear;

% Input the E.coli core model
model=readCbModel('iECD_1391');

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
model = changeRxnBounds(model,'EX_o2(e)',-2,'l');
% model = changeRxnBounds(model,'EX_biliverdin(e)',0.00322,'l'); % Estimated Chen values
model = changeRxnBounds(model,'EX_biliverdin(e)',0.05,'l'); % Large production
model = changeObjective(model,'Ec_biomass_iJO1366_core_53p95M');

for i=1:10
     model = changeRxnBounds(model,'EX_o2(e)',-2*i,'l');
     FBAsolution = optimizeCbModel(model,'max',0,0);
%      secretionTable(i) = printFluxVector(model,FBAsolution.x, true, true);
     printFluxVector(model,FBAsolution.x, true, true);
end
