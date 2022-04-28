% EssentialGenes.m
clear; clc;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

tol = 1e-6; % Growth rate lower limit
grRatio = singleGeneDeletion(model);
grRatio(isnan(grRatio))=0;
EssentialGenes = model.genes(grRatio<tol)

[results ListResults] = findRxnsFromGenes(model, EssentialGenes);

