% DynamicGrowth_Aerobic_JO1366.m
clear;

model = readCbModel('iJO1366.mat');
model = changeRxnBounds(model, {'EX_glc(e)','EX_o2(e)'},[-10 -5], 'l'); 
model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');

% Set-up variables for dynamicFBA
substrateRxns = {'EX_glc(e)'};
initConcentrations = [10]; 
initBiomass = .01;
timeStep = .25; nSteps = 100;
plotRxns = {'EX_ac(e)','EX_acald(e)','EX_etoh(e)','EX_for(e)','EX_glc(e)','EX_lac-L(e)','EX_succ(e)'};

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