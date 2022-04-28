% Test_modelReg.m
clear; clc;

model = readCbModel('modelReg.mat');

% Initial Conditions
model = changeRxnBounds(model, 'EX_glc(e)',-10, 'l');
model = changeRxnBounds(model, 'EX_fum(e)',-10, 'l');
model = changeRxnBounds(model, 'EX_o2(e)',-30, 'l');

% Set optimization objective
model = changeObjective(model,'Biomass_Ecoli_core_w_GAM');

% FBA analysis using regulated constraints
[FBAsols,DRgenes,constrainedRxns,cycleStart,states]= optimizeRegModel(model);

% Create a statevector showing the states created in teh regulatory process
stateVector = [model.regulatoryGenes;model.regulatoryInputs1;model.regulatoryInputs2];

% Print the regulatory states
printLabeledData(stateVector,[states{1,1},states{1,2},states{1,3},states{1,4},states{1,5},states{1,6},states{1,7}])

% Print the regulated optimized FBA fluxes
printFluxVector(model, FBAsols{1,1}.x, true)

% Draw the flux values on the map "target.svg"
map=readCbMap('ecoli_Textbook_ExportMap');
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;
drawFlux(map, model, FBAsols{1,1}.x, options);