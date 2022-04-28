% DynamicBiliverdinProduction_Nutrient_Uptake_BL21_Plasmid.m
clear;
changeCobraSolver('glpk','all'); % LP solver set to glpk 

model=readCbModel('iECD_1391.xml');
model = removeRxns(model,'Ec_biomass_iJO1366_WT_53p95M')

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
rxnID=findRxnIDs(model,'EX_biliverdin(e)');
model.subSystems(rxnID) = {'Exchange'};

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

model = changeRxnBounds(model, {'EX_glc(e)','EX_o2(e)'},[-10 -20], 'l'); % Aerobic
% model = changeRxnBounds(model,'EX_biliverdin(e)',0,'l'); % Biliverdin Exchange Reaction
model = changeRxnBounds(model,'HEMEOX',0.1,'l'); 
model = changeObjective(model,'Ec_biomass_iJO1366_core_53p95M');

% Set-up  for one variables for dynamicFBA
% NOTE- substrate rxns and plot rxns need to be in the order that they appear in the model

% K-12 Concentrations
substrateRxns = {'EX_cl(e)','EX_cu2(e)','EX_fe3(e)','EX_glc(e)','EX_k(e)','EX_mg2(e)','EX_mn2(e)','EX_mobd(e)',...
    'EX_na1(e)','EX_nh4(e)','EX_pi(e)','EX_so4(e)','EX_zn2(e)'};
initConcentrations = [1.750214773,0.008010093,0.087865335,138.7655417,49.74894838,2.028562155,0.101058487,0.019031286,...
    1.295164976,75.71742258,70.08147994,2.11641021,0.017388304];
plotRxns = {'EX_ala-L(e)','EX_arg-L(e)','EX_asn-L(e)','EX_asp-L(e)','EX_cl(e)','EX_cu2(e)','EX_cys-L(e)','EX_fe3(e)',...
    'EX_glc(e)','EX_gln-L(e)','EX_glu-L(e)','EX_gly(e)','EX_his-L(e)','EX_ile-L(e)','EX_k(e)','EX_leu-L(e)',...
    'EX_lys-L(e)','EX_met-L(e)','EX_mg2(e)','EX_mn2(e)','EX_mobd(e)','EX_na1(e)','EX_nh4(e)','EX_phe-L(e)',...
    'EX_pi(e)','EX_pro-L(e)','EX_ser-L(e)','EX_so4(e)','EX_thm(e)','EX_thr-L(e)','EX_tyr-L(e)','EX_val-L(e)','EX_zn2(e)'};

initBiomass = 0.01; % g/L
timeStep = .25; nSteps = 500;
plotRxns = {'EX_glc(e)','EX_k(e)','EX_nh4(e)','EX_pi(e)','EX_so4(e)'}; 

figure(1)
[concentrationMatrix,excRxnNames,timeVec,biomassVec] = ...
    dynamicFBA(model,substrateRxns,initConcentrations, initBiomass, timeStep, nSteps, plotRxns);

% Plot labels
subplot(1,2,1); title('Biomass Concentration'); xlabel('Time(h)'); ylabel('Concentration (g/L)');
subplot(1,2,2); title('Extracellular Concentrations'); xlabel('Time(h)'); ylabel('Concentrations (mmol/L)');

cMatrix = full(concentrationMatrix);
