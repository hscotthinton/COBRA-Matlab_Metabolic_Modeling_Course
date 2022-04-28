% succinateOptimalSolutionsKnockouts.m
clear; clc;

model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',0,'l');
model = changeRxnBounds(model,'EX_o2(e)',-40,'l');
model = changeRxnBounds(model,'EX_succ(e)',-20,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Delate reactions to select desired alternate flux vector
model = changeRxnBounds(model,'ME1',0,'b'); % Knockout the ME1 reaction
model = changeRxnBounds(model,'PYK',0,'b'); % Knockout the PYK reaction

% List optimal solutions
% solverOK = changeCobraSolver('glpk','MILP'); % Won't work with Gurobi
[solutions] = enumerateOptimalSolutions(model); 
