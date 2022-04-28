% Biliverdin_Carbon_Source_BL21_Plasmid.m
clear;



% Input the E.coli core model
model=readCbModel('iECD_1391');
model = removeRxns(model,'Ec_biomass_iJO1366_WT_53p95M'); % Remove WT biomass function

% Add plasmid nucleotide precursors
model=addReaction(model,'PLASMID','3276.86 dgtp[c] + 3276.86 dctp[c] + 3199.14 datp[c] + 3199.14 dttp[c] -> plasmid[c]');

% Add demand reaction for plasmid DNA and set flux rate
model = addDemandReaction(model,'plasmid[c]');
model = changeRxnBounds(model,'DM_plasmid[c]',0.63e-9,'l');

% Add beta-lactamase
model=addReaction(model,'b-lactamase','28 ala-L[c] + 3 cys-L[c] + 16 asp-L[c] + 20 glu-L[c] + 9 phe-L[c] + 21 gly[c] + 7 his-L[c] + 17 ile-L[c] + 11 lys-L[c] + 33 leu-L[c] + 10 met-L[c] + 8 asn-L[c] + 14 pro-L[c] + 9 gln-L[c] + 19 arg-L[c] + 17 ser-L[c] + 20 thr-L[c] + 16 val-L[c] + 4 trp-L[c] + 4 tyr-L[c] + 1,231.52 atp[c] -> b-lactamase[c] + 1,231.52 adp[c] + 1,231.52 pi[c]');

% Add demand reaction for beta-lactamase
model = addDemandReaction(model,'b-lactamase[c]');
model = changeRxnBounds(model,'DM_b-lactamase[c]',0.000569,'l');

% Add heme oxygenase enzyme
model=addReaction(model,'HOprotein','35 ala-L[c] + 3 cys-L[c] + 14 asp-L[c] + 36 glu-L[c] + 17 phe-L[c] + 19 gly[c] + 6 his-L[c] + 6 ile-L[c] + 26 lys-L[c] + 32 leu-L[c] + 13 met-L[c] + 14 asn-L[c] + 10 pro-L[c] + 16 gln-L[c] + 13 arg-L[c] + 12 ser-L[c] + 14 thr-L[c] + 13 val-L[c] + 2 trp-L[c] + 14 tyr-L[c] + 1,356.39 atp[c] -> HO[c] + 1,356.39 adp[c] + 1,356.39 pi[c]');

% Add demand reaction for heme oxygenase enzyme
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
model = changeRxnBounds(model,'Ec_biomass_iJO1366_WT_53p95M',0,'b');

% Set key variables
model = changeObjective(model,'EX_biliverdin(e)');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
number = 1;
disp(' ');

% Optimization objective to glucose
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
% [carbonSources(number)] = {'Glucose', FBAsolution.f}
disp({'Glucose', FBAsolution.f});

% Optimization objective to glycerol
model = changeRxnBounds(model,'EX_glc(e)',-0,'l');
model = changeRxnBounds(model,'EX_glyc(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Glycerol', FBAsolution.f});

% Optimization objective to lactose
model = changeRxnBounds(model,'EX_glyc(e)',-0,'l');
model = changeRxnBounds(model,'EX_lcts(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Lactose', FBAsolution.f});

% Optimization objective to sucrose
model = changeRxnBounds(model,'EX_lcts(e)',-0,'l');
model = changeRxnBounds(model,'EX_sucr(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Sucrose', FBAsolution.f});

% Optimization objective to fructose
model = changeRxnBounds(model,'EX_sucr(e)',-0,'l');
model = changeRxnBounds(model,'EX_fru(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Fructose', FBAsolution.f});

% Optimization objective to fructose
model = changeRxnBounds(model,'EX_fru(e)',-0,'l');
model = changeRxnBounds(model,'EX_malt(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Maltose', FBAsolution.f});

% Optimization objective to D-mannitol
model = changeRxnBounds(model,'EX_malt(e)',-0,'l');
model = changeRxnBounds(model,'EX_mnl(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'D-mannitol', FBAsolution.f});

% Optimization objective to D-sorbitol
model = changeRxnBounds(model,'EX_mnl(e)',-0,'l');
model = changeRxnBounds(model,'EX_sbt-D(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'D-sorbitol', FBAsolution.f});

% Optimization objective to succinate
model = changeRxnBounds(model,'EX_sbt-D(e)',-0,'l');
model = changeRxnBounds(model,'EX_succ(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Succinate', FBAsolution.f});

% Optimization objective to malate
model = changeRxnBounds(model,'EX_succ(e)',-0,'l');
model = changeRxnBounds(model,'EX_mal-D(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Malate', FBAsolution.f});

% Optimization objective to citrate
model = changeRxnBounds(model,'EX_mal-D(e)',-0,'l');
model = changeRxnBounds(model,'EX_cit(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Citrate', FBAsolution.f});

% Optimization objective to L-glutamate
model = changeRxnBounds(model,'EX_cit(e)',-0,'l');
model = changeRxnBounds(model,'EX_glu-L(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'L-glutamate', FBAsolution.f});

% Optimization objective to L-glutamate + glucose
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_glu-L(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Glucose + L-glutamate', FBAsolution.f});

% Optimization objective to Lactose + Glycerol
model = changeRxnBounds(model,'EX_glc(e)',-0,'l');
model = changeRxnBounds(model,'EX_glu-L(e)',-0,'l');
model = changeRxnBounds(model,'EX_lcts(e)',-10,'l');
model = changeRxnBounds(model,'EX_glyc(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Lactose + Glycerol', FBAsolution.f});

% Optimization objective to Mannose
model = changeRxnBounds(model,'EX_glc(e)',-0,'l');
model = changeRxnBounds(model,'EX_lcts(e)',-0,'l');
model = changeRxnBounds(model,'EX_glyc(e)',-0,'l');
model = changeRxnBounds(model,'EX_man(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Mannose', FBAsolution.f});

% Optimization objective to Rhamnose
model = changeRxnBounds(model,'EX_glc(e)',-0,'l');
model = changeRxnBounds(model,'EX_man(e)',-0,'l');
model = changeRxnBounds(model,'EX_rmn(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Rhamnose', FBAsolution.f});

% Optimization objective to Ribose
model = changeRxnBounds(model,'EX_glc(e)',-0,'l');
model = changeRxnBounds(model,'EX_rmn(e)',-0,'l');
model = changeRxnBounds(model,'EX_rib-D(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'Ribose', FBAsolution.f});

% Optimization objective to xylose
model = changeRxnBounds(model,'EX_glc(e)',-0,'l');
model = changeRxnBounds(model,'EX_rib-D(e)',-0,'l');
model = changeRxnBounds(model,'EX_xyl-D(e)',-10,'l');
FBAsolution = optimizeCbModel(model,'max');
disp({'xylose', FBAsolution.f});