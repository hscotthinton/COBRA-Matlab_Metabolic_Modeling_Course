% Sampling_Histogram.m
clear;

% Load the E.coli core model and set constraints
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Sample model
[sampleStruct,mixedFrac] = gpSampler(model,5000,[],120);
% [sampleStruct,mixedFrac] = gpSampler(model,5000,[],[],1000);

% Determine the minimum and maximum possible fluxes so the sampling results 
% can be plotted for the first nine reactions in the model

[minFlux,maxFlux] = fluxVariability(model,0);
figure(1);
for i = 1 : 95
    subplot(8,12,i)
    hist(sampleStruct.points(i,:),50);
    hold on
    plot([minFlux(i) maxFlux(i)], [0 1],'*r');
    title(model.rxns{i});
end

