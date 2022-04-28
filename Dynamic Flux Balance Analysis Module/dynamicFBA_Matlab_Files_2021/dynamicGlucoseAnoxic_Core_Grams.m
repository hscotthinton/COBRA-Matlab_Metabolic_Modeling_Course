% dynamicGlucoseAnoxic_Core_Grams.m
clear; 

model = readCbModel('ecoli_core_model.mat');

model = changeRxnBounds(model, 'EX_glc(e)',-10, 'l');
model = changeRxnBounds(model, 'EX_o2(e)',-10, 'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

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
subplot(1,2,2); title('Substrate Concentrations'); xlabel('Time steps'); ylabel('Concentrations (mmol/L)');

CM = full(concentrationMatrix);

% Plot formate production in mg/L
CM = full(concentrationMatrix);
% % [a,floc] = ismember('EX_for(e)',excRxnNames);
% % [a,aloc] = ismember('EX_ac(e)',excRxnNames);
[a,loc] = ismember({'EX_for(e)','EX_ac(e)'},plotRxns);

formateMW = 0.045; % formate  molecular weight (g/mmol)
acetateMW = 0.236; % acetate  molecular weight (g/mmol)

figure(2)
clf;
subplot(1,2,1); 
plot(timeVec,biomassVec);
title('Biomass Concentration'); xlabel('Time (h)'); ylabel('Concentration (g/L)');
axis tight
subplot(1,2,2); 
% plot(timeVec,formateMW*concentrationMatrix(loc,:));
plot(timeVec,formateMW*concentrationMatrix(loc,:));
title('Substrate Concentrations'); xlabel('Time (h)'); ylabel('Concentration (g/L)');
axis tight
legend(strrep(excRxnNames(loc),'EX_',''));