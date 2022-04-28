% DynamicBiliverdinProduction_Cont_Plasmid_BL21.m
clear;

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
model = changeObjective(model,'Ec_biomass_iJO1366_core_53p95M');

model = changeRxnBounds(model,'HEMEOX',0.8,'l'); 
model = changeRxnBounds(model,'EX_biliverdin(e)',0,'l'); 

% Set-up  for one variables for dynamicFBA
substrateRxns = {'EX_glc(e)'};
initConcentrations = [50];
initBiomass = .01;
timeStep = .25; nSteps = 32; % nSteps = 32 for 8 hours;
plotRxns = {'EX_biliverdin(e)'};
% plotRxns = {'EX_ac(e)','EX_etoh(e)','EX_for(e)','EX_glyclt(e)','DM_HO[c]','EX_biliverdin(e)'};
exclUptakeRxns = {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)','EX_nh4(e)','EX_pi(e)','EX_so4(e)','EX_k(e)'};

[concentrationMatrix,excRxnNames,timeVec,biomassVec] = ...
    dynamicFBA(model,substrateRxns,initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns );

% Plot labels
subplot(1,2,1); title('Biomass Concentration'); xlabel('Time (h)'); ylabel('Concentration (g/L)');
subplot(1,2,2); title('Extracellular Concentrations'); xlabel('Time (h)'); ylabel('Concentration (mmol/L)');

CM = full(concentrationMatrix);

% Plot biliverdin production in mg/L
CM = full(concentrationMatrix);
[a,loc] = ismember('EX_biliverdin(e)',excRxnNames)
biliverdinMW = 0.582646; % Biliverdin molecular weight (mmol/g)

figure(2)
clf;
subplot(1,2,1); 
plot(timeVec,biomassVec);
subplot(1,2,1); title('Biomass Concentration'); xlabel('Time (h)'); ylabel('Concentration (g/L)');
axis tight
subplot(1,2,2); 
plot(timeVec,biliverdinMW*concentrationMatrix(loc,:));
% title('Substrate Concentrations'); xlabel('Time (h)'); ylabel('Concentration (mg/L)');
subplot(1,2,2); title('Extracellular Concentrations'); xlabel('Time (h)'); ylabel('Concentration (g/L)');
axis tight
legend(strrep(excRxnNames(loc),'EX_',''));
