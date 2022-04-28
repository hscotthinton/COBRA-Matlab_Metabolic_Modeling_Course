% findingOptimalSolutionsEthanol.m
clear;

model = readCbModel('ecoli_core_model.mat');

model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

model = changeObjective(model,'EX_etoh(e)');

changeCobraSolver('glpk','MILP'); % Won't work with Gurobi

% List optimal solutions
[solutions] = enumerateOptimalSolutions(model); 
