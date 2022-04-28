% IdentifyingCorrelSets.m
clear;

% Input the E.coli core model and set constraints
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

[sampleStruct,mixedFrac] = gpSampler(model,5000,[],120);
[setsSorted,setNoSorted,setSize] = identifyCorrelSets(model,sampleStruct.points);
setNames = [];
setNumbers = [];
disp('Correlated Reations Sets')
for i = 1 : length(setsSorted)
    setNames = [setNames; setsSorted{i}.names];
    setNumbers = [setNumbers;i*ones(length(setsSorted{i}.names),1)];
    disp([i,setsSorted{i}.names'])
end