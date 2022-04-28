% Sampling_Succ_Optimal.m
clear; 

model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',0,'l');
model = changeRxnBounds(model,'EX_o2(e)',-40,'l');
model = changeRxnBounds(model,'EX_succ(e)',-20,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Find optimal value objective function (growth rate)
FBAsolution = optimizeCbModel(model,'max',0,0);
printFluxVector(model, FBAsolution.x, true);

% Force output to be close to optimal value
% model = changeRxnBounds(model,'Biomass_Ecoli_core_N(w/GAM)_Nmet2',FBAsolution.f,'b');
model = changeRxnBounds(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2',0.9*FBAsolution.f,'l');

% Sample model
[sampleStruct,mixedFrac] = gpSampler(model,5000,[],120);

% Plot histograms for selected reactions
rxnList = {'FORt2', 'FORti', 'MDH', 'ME1', 'ME2', 'NADTRHD', 'PPCK', 'PYK', 'EX_succ(e)', 'Biomass_Ecoli_core_N(w/GAM)-Nmet2'};

% Include flux variability analysis mins and maxs on the histogram
figure(1);
[minFlux,maxFlux]=fluxVariability(model,100,'max',model.rxns,false,false);
for i = 1 : 10
    subplot(2,5,i)
    rxnID = findRxnIDs(model,rxnList(i)); 
    hist(sampleStruct.points(rxnID,:),50);
    hold on
    plot([minFlux(rxnID) maxFlux(rxnID)], [0 1],'*r');
    title(rxnList(i));
end