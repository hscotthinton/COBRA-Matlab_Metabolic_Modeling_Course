% findingOptimalSolutionsSucc.m
clear; clc;

model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',0,'l');
model = changeRxnBounds(model,'EX_o2(e)',-40,'l');
model = changeRxnBounds(model,'EX_succ(e)',-20,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% List optimal solutions
% changeCobraSolver('glpk','all'); % Won't work with Gurobi
[solutions] = enumerateOptimalSolutions(model); 
