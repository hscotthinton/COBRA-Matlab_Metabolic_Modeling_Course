% Catabolite_Repression_iMC1010v2.m
clear; 

model = readCbModel('iMC1010v2.mat'); % Reconstructed model found in lesson Matlab files

% Begin analysis
model = changeRxnBounds(model, 'EX_glc(e)',-10, 'l');
model = changeRxnBounds(model, 'EX_lcts(e)',-0, 'l');
model = changeRxnBounds(model, 'EX_glyc(e)',-10, 'l');
model = changeRxnBounds(model, 'EX_succ(e)',-0, 'l');
model = changeRxnBounds(model, 'EX_fum(e)',-0, 'l');
model = changeRxnBounds(model, 'EX_o2(e)',-5, 'l');

[FBAsols,DRgenes,constrainedRxns,cycleStart,states]= optimizeRegModel(model);

printFluxVector(model, FBAsols{1,1}.x, true)

substrateRxns = {'EX_fum(e)','EX_glc(e)','EX_glyc(e)','EX_lcts(e)','EX_succ(e)'};
initConcentrations = [0 10 15 0 0];

initBiomass = .035;
timeStep = .25; 
nSteps = 100;
plotRxns = {'EX_ac(e)','EX_co2(e)','EX_etoh(e)','EX_for(e)','EX_fru(e)','EX_fum(e)','EX_glc(e)',...
    'EX_gln_L(e)','EX_glu_L(e)','EX_glyc(e)','EX_lac_D(e)','EX_lcts(e)','EX_mal_L(e)','EX_pro_L(e)','EX_pyr(e)','EX_succ(e)'};

[concentrationMatrix,excRxnNames,timeVec,biomassVec,drGenes,constrainedRxns,states] = ...
    dynamicRFBA(model,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,plotRxns);

cMatrix = full(concentrationMatrix);

% Plot labels
subplot(1,2,1); title('Biomass Concentration'); xlabel('Time steps'); ylabel('Concentration (g/L)');
subplot(1,2,2); title('Extracellular Concentrations'); xlabel('Time steps'); ylabel('Concentrations (mmol/L)');