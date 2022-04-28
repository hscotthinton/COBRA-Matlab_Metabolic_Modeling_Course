% Biliverdin_Robustness_HEMEOX_BL21.m
clear;

% Input the E.coli core model
model=readCbModel('iECD_1391');

model = removeRxns(model,'Ec_biomass_iJO1366_WT_53p95M'); % remove wild type biomass function

% Add heme oxygenase enzyme and Demand Function
model=addReaction(model,'HOprotein','35 ala-L[c] + 3 cys-L[c] + 14 asp-L[c] + 36 glu-L[c] + 17 phe-L[c] + 19 gly[c] + 6 his-L[c] + 6 ile-L[c] + 26 lys-L[c] + 32 leu-L[c] + 13 met-L[c] + 14 asn-L[c] + 10 pro-L[c] + 16 gln-L[c] + 13 arg-L[c] + 12 ser-L[c] + 14 thr-L[c] + 13 val-L[c] + 2 trp-L[c] + 14 tyr-L[c] + 1,356.39 atp[c] -> HO[c] + 1,356.39 adp[c] + 1,356.39 pi[c]');
model = addDemandReaction(model,'HO[c]');
model = changeRxnBounds(model,'DM_HO[c]',0.0102,'l');

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

% Set optimization objective
model = changeObjective(model,'Ec_biomass_iJO1366_core_53p95M');

% Robustenss Analysis
[controlFlux, objFlux] = robustnessAnalysis(model,'EX_biliverdin(e)',100);