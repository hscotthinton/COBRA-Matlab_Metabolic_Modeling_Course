% Glucose_Low_Aerobic_iMC1010v2.m
clear; 


model = readCbModel('iMC1010v2.mat'); % Reconstructed model found in lesson Matlab files

model = changeRxnBounds(model, 'EX_glc(e)',-10, 'l');
model = changeRxnBounds(model, 'EX_o2(e)',-5, 'l');

[FBAsols,DRgenes,constrainedRxns,cycleStart,states]= optimizeRegModel(model);

substrateRxns = {'EX_glc(e)'};
initConcentrations = [20]; 
initBiomass = .035;
timeStep = .25; 
nSteps = 100;
plotRxns = {'EX_ac(e)','EX_etoh(e)','EX_for(e)','EX_fru(e)','EX_glc(e)',...
    'EX_gln_L(e)','EX_glu_L(e)','EX_lac_D(e)','EX_mal_L(e)','EX_succ(e)'};

[concentrationMatrix,excRxnNames,timeVec,biomassVec,drGenes,constrainedRxns,states] = ...
    dynamicRFBA(model,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,plotRxns);

% Plot labels
subplot(1,2,1); title('Biomass Concentration'); xlabel('Time steps'); ylabel('Concentration (g/L)');
subplot(1,2,2); title('Substrate Concentrations'); xlabel('Time steps'); ylabel('Concentrations (mmol/L)');