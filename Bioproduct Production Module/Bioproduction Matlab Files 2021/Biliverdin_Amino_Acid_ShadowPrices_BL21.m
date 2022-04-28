% Biliverdin_Amino_Acid_ShadowPrices_BL21.m
clear;
% changeCobraSolver('glpk','all'); % LP solver set to glpk 
changeCobraSolverParams('LP', 'primalOnly', false); % Required for shadow prices and reduced costs

% Input the E.coli core model
model=readCbModel('iECD_1391');
model = removeRxns(model,'Ec_biomass_iJO1366_WT_53p95M');

% Add plasmid nucleotide precursors
model=addReaction(model,'PLASMID','3276.86 dgtp[c] + 3276.86 dctp[c] + 3199.14 datp[c] + 3199.14 dttp[c] -> plasmid[c]');
model = changeRxnBounds(model,'PLASMID',0.63e-9,'b');

% Add demand reaction for plasmid DNA and set flux rate
model = addDemandReaction(model,'plasmid[c]'); %DM_plasmid[c]

% Add beta-lactamase
model=addReaction(model,'b-lactamase','28 ala-L[c] + 3 cys-L[c] + 16 asp-L[c] + 20 glu-L[c] + 9 phe-L[c] + 21 gly[c] + 7 his-L[c] + 17 ile-L[c] + 11 lys-L[c] + 33 leu-L[c] + 10 met-L[c] + 8 asn-L[c] + 14 pro-L[c] + 9 gln-L[c] + 19 arg-L[c] + 17 ser-L[c] + 20 thr-L[c] + 16 val-L[c] + 4 trp-L[c] + 4 tyr-L[c] + 1,231.52 atp[c] -> b-lactamase[c] + 1,231.52 adp[c] + 1,231.52 pi[c]');
model = changeRxnBounds(model,'b-lactamase',0.000569,'b');

% Add demand reaction for beta-lactamase
model = addDemandReaction(model,'b-lactamase[c]'); % DM_B-lactamase[c]

% Add heme oxygenase enzyme
model=addReaction(model,'HOprotein','35 ala-L[c] + 3 cys-L[c] + 14 asp-L[c] + 36 glu-L[c] + 17 phe-L[c] + 19 gly[c] + 6 his-L[c] + 6 ile-L[c] + 26 lys-L[c] + 32 leu-L[c] + 13 met-L[c] + 14 asn-L[c] + 10 pro-L[c] + 16 gln-L[c] + 13 arg-L[c] + 12 ser-L[c] + 14 thr-L[c] + 13 val-L[c] + 2 trp-L[c] + 14 tyr-L[c] + 1,356.39 atp[c] -> HO[c] + 1,356.39 adp[c] + 1,356.39 pi[c]');
model = changeRxnBounds(model,'HOprotein',0.0102,'b');

% Add demand reaction for heme oxygenase enzyme
model = addDemandReaction(model,'HO[c]'); % DM_HO[c]

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

% Optimization objective to maximize bilverdin production
model = changeObjective(model,'EX_biliverdin(e)');
% model = changeObjective(model,'Ec_biomass_iJO1366_core_53p95M');
FBAsolution = optimizeCbModel(model,'max');

'Flux values'
% printFluxVector(model,FBAsolution.x, true)

'Shadow prices'
% printShadowPriceVector(model, FBAsolution.y)

aminoAcids = {'ala-L[e]','arg-L[e]','asn-L[e]','asp-L[e]','cys-L[e]','gln-L[e]','glu-L[e]',...
    'gly[e]','h2o[e]','his-L[e]','ile-L[e]','leu-L[e]','lys-L[e]','met-L[e]','phe-L[e]',...
    'pro-L[e]','ser-L[e]','thr-L[e]','trp-L[e]','tyr-L[e]','val-L[e]','pi[e]','cbl1[e]'};

[tmp,preNumber] = ismember(aminoAcids,model.mets);
aminoAcidsSP = FBAsolution.y(preNumber);
metIncrease = 1./aminoAcidsSP;

printLabeledData(aminoAcids',aminoAcidsSP)

% 'Number of metabolites required to increase biliverdin prduction by one'
% printLabeledData(precursors',metIncrease)