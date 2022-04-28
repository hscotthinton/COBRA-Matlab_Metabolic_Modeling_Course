% EthanolProduction_multiProductionEnvelope_iaf1260.m
clear; 

model = readCbModel('iAF1260.mat');

model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');
model = changeObjective(model,'Ec_biomass_iAF1260_core_59p81M');

% deletions = {}; 
deletions = {'GLUDy','TPI','PFL'};  

biomassRxn = {'Ec_biomass_iAF1260_core_59p81M'};

%Show only growth coupled metabolites
figure(1)
[biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,20,false);

% %Show all secreted metabolites
% figure(2)
% [biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,20,true);