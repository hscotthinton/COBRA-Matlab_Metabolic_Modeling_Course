% EthanolProduction_simpleOptKnock.m
clear; clc;

model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% deletions = {};
% deletions = {'PFL'}; 
% deletions = {'PTAr','PYK'}; 
% deletions = {'ACKr','GLUDy','PYK'}; 
% deletions = {'GND','ME2','PTAr','PYK'}; 
deletions = {'ACKr','GLUDy','ME2','PGL','PYK'}; 

biomassRxn = {'Biomass_Ecoli_core_N(w/GAM)-Nmet2'};
minGrowth = 0.05;
geneDelFlag = false;
doubleDelFlag = false;

[wtRes,delRes] = simpleOptKnock(model,'EX_etoh(e)',deletions,false,minGrowth,doubleDelFlag) 