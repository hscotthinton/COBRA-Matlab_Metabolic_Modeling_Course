% EthanolProduction_doubleProductionEnvelope.m
clear; clc;

model = readCbModel('ecoli_core_model.mat');

model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

deletions = {};
% deletions = {'PFL'}; 
% deletions = {'PTAr','PYK'}; 
% deletions = {'ACKr','GLUDy','PYK'}; 
% deletions = {'GND','ME2','PTAr','PYK'}; 
% deletions = {'ACKr','GLUDy','ME2','PGL','PYK'}; 

biomassRxn = {'Biomass_Ecoli_core_N(w/GAM)-Nmet2'};

[x1,x2,y] = doubleProductionEnvelope(model,deletions,'EX_etoh(e)','EX_for(e)',biomassRxn,false,20);
