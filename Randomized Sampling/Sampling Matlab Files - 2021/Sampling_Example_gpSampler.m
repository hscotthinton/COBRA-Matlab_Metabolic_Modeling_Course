% Sampling_Example_gpSampler.m
clear;

% Input the E.coli core model and set constraints
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
% model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

model_aerobic = model;
sampleStruct_aerobic = gpSampler(model_aerobic,2000,[],120);
% Simulation time is ~120 s with the Gurobi LP solver.
model_anaerobic = changeRxnBounds(model_aerobic,'EX_o2(e)',0,'l');
sampleStruct_anaerobic = gpSampler(model_anaerobic,2000,[],120);
% Simulation time is ~120 s with the Gurobi LP solver. Sampling results will be returned in the two structures sampleStruct_aerobic and sampleStruct_anaerobic within the field points.
% Visualize sampling results for a set of reactions. 
rxnList = {'PGI', 'PFK', 'FBP', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'};
plotSampleHist(rxnList, {sampleStruct_aerobic.points, sampleStruct_anaerobic.points }, {model_aerobic, model_anaerobic},[],[2,5]);