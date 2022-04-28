% DynamicGrowth_Aerobic_JO1366.m
clear;

model = readCbModel('iJO1366.mat');
model = changeRxnBounds(model, {'EX_glc__D_e','EX_o2_e'},[-10 -10], 'l'); 
model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');

% Set-up variables for dynamicFBA
substrateRxns = {'EX_glc__D_e'};
initConcentrations = [10]; 
initBiomass = .01;
timeStep = .25; nSteps = 100;
plotRxns = {'EX_ac_e','EX_acald_e','EX_etoh_e','EX_for_e','EX_glc__D_e','EX_lac__L_e','EX_succ_e'};

[concentrationMatrix,excRxnNames,timeVec,biomassVec] = ...
    dynamicFBA(model,substrateRxns,initConcentrations, initBiomass, timeStep, nSteps, plotRxns);

% Plot labels
subplot(1,2,1);
title('Biomass Concentration');
xlabel('Time steps');
ylabel('Concentration (g/L)');
subplot(1,2,2);
title('Extracellular Concentrations');
xlabel('Time steps');
ylabel('Metabolite Concentrations (mmol/L)');