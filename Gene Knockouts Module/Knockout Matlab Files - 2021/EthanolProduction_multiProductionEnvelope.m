% EthanolProduction_multiProductionEnvelope.m
clear; 

model = readCbModel('ecoli_core_model.mat');

model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

% deletions = {}; 
% deletions = {'PFL'}; 
% deletions = {'PTAr','G6PDH2r'}; 
% deletions = {'ACKr','GLUDy','G6PDH2r'}; 
% deletions = {'ACKr','GLUDy','PYK'}; 
deletions = {'ACKr','GLUDy','ME2','PGL','PYK'}; 

biomassRxn = {'Biomass_Ecoli_core_N(w/GAM)-Nmet2'};

%Show only growth coupled metabolites
figure(1)
[biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,20,false);

%Show all secreted metabolites
figure(2)
[biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,20,true);