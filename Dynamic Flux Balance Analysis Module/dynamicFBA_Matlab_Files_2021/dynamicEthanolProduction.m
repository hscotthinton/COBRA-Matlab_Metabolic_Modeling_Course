% dynamicEthanolProduction.m
clear;

% Load the E.coli core model
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% 3 Knockout reactions ( Growth-rate > 0.05 )
model = changeRxnBounds(model,'PTAr',0,'b');
model = changeRxnBounds(model,'GLUDy',0,'b');
model = changeRxnBounds(model,'PYK',0,'b');

% Set-up variables for dynamicFBA
substrateRxns = {'EX_glc(e)'};
initConcentrations = [10]; initBiomass = .035;
timeStep = .25; nSteps = 100;
plotRxns = {'EX_ac(e)','EX_etoh(e)','EX_for(e)','EX_fru(e)','EX_fum(e)','EX_glc(e)','EX_gln_L(e)',...
    'EX_glu_L(e)','EX_lac_D(e)','EX_mal_L(e)','EX_succ(e)','EX_acald(e)'};

 [concentrationMatrix,excRxnNames,timeVec,biomassVec] = ...
     dynamicFBA(model,substrateRxns,initConcentrations, initBiomass, timeStep, nSteps, plotRxns);

 % Plot labels
subplot(1,2,1); title('Biomass Concentration'); xlabel('Time steps'); ylabel('Concentration (g/L)');
subplot(1,2,2); title('Extracellular Concentrations'); xlabel('Time steps'); ylabel('Metabolite Concentrations (mmol/L)');

% Plot ethanol production in g/L
CM = full(concentrationMatrix);
[a,loc] = ismember({'EX_etoh(e)'},excRxnNames);
ethanolMW = 0.04607; % ethanol molecular weight (g/mmol)

figure(2)
clf;
subplot(1,2,1); 
plot(timeVec,biomassVec);
title('Biomass Concentration'); xlabel('Time (h)'); ylabel('Concentration (g/L)');
axis tight
subplot(1,2,2); 
plot(timeVec,ethanolMW*concentrationMatrix(loc,:));
title('Extracellular Concentrations'); xlabel('Time (h)'); ylabel('Metabolite Concentration (g/L)');
axis tight
legend(strrep(excRxnNames(loc),'EX_',''));