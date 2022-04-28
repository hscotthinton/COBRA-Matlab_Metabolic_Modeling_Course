% GlucoseFructose_Aerobic.m
clear; 

model = readCbModel('modelReg.mat');
model = changeRxnBounds(model, 'EX_glc(e)',-10, 'l');
model = changeRxnBounds(model, 'EX_fru(e)',-10, 'l'); % Normally -10
model = changeRxnBounds(model, 'EX_o2(e)',-5, 'l'); % Normally -30

% Set optimization objective
model = changeObjective(model,'Biomass_Ecoli_core_w_GAM');

% [FBAsols,DRgenes,constrainedRxns,cycleStart,states]= optimizeRegModel(modelReg);

substrateRxns = {'EX_fru(e)','EX_glc(e)'};
initConcentrations = [5 10]; % Normally [5 10]
initBiomass = .035;
timeStep = .25; 
nSteps = 100;
plotRxns = {'EX_ac(e)','EX_etoh(e)','EX_for(e)','EX_fru(e)','EX_glc(e)','EX_gln_L(e)',...
    'EX_glu_L(e)','EX_lac_D(e)','EX_mal_L(e)','EX_succ(e)'};

% plotRxns = {'EX_ac(e)','EX_etoh(e)','EX_for(e)','EX_fru(e)','EX_glc(e)','EX_gln_L(e)',...
%     'EX_glu_L(e)','EX_lac_D(e)','EX_mal_L(e)','EX_succ(e)','EX_acald(e)'};

[concentrationMatrix,excRxnNames,timeVec,biomassVec,drGenes,constrainedRxns,states] = ...
    dynamicRFBA(model,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,plotRxns);

cMatrix = full(concentrationMatrix);

% Plot labels
subplot(1,2,1); title('Biomass Concentration'); xlabel('Time steps'); ylabel('Concentration (g/L)');
subplot(1,2,2); title('Substrate Concentrations'); xlabel('Time steps'); ylabel('Concentrations (mmol/L)');