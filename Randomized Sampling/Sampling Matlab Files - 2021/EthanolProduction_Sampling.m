% EthanolProduction_Sampling.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set carbon source and oxygen uptake rates for wild type model
model = changeRxnBounds(model,'EX_glc(e)',-5,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');
FBAsolution = optimizeCbModel(model,'max',0,0);
model_WT = model;

% Knockout reactions for mutant model
model = changeRxnBounds(model,'NADH16',0,'b');
model = changeRxnBounds(model,'PTAr',0,'b');
model = changeRxnBounds(model,'TKT2',0,'b');
Mutantsolution = optimizeCbModel(model,'max',0,0);
model_Mutant = model;

sampleStruct_WT = gpSampler(model_WT,2000,[],120);
% Simulation time is ~120 s with the Gurobi LP solver.

sampleStruct_Mutant = gpSampler(model_Mutant,2000,[],120);
% Simulation time is ~120 s with the Gurobi LP solver. Sampling results will be returned in the two
% structures sampleStruct_WT and sampleStruct_Mutant within the field points.
% Visualize sampling results for a set of reactions. 
rxnList = {'EX_glc(e)', 'EX_o2(e)', 'EX_etoh(e)', 'EX_lac_D(e)', 'EX_for(e)', 'ATPS4r', 'CYTBD', 'GND', 'ICDHyr', 'Biomass_Ecoli_core_N(w/GAM)-Nmet2'};
plotSampleHist(rxnList, {sampleStruct_WT.points, sampleStruct_Mutant.points }, {model_WT, model_Mutant},[],[2,5]);