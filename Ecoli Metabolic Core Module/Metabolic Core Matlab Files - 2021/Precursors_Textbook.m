% ATP_Textbook.m

clear;

model = readCbModel('ecoli_core_model.mat');

model = changeRxnBounds(model,'EX_glc(e)',-1,'l');
model = changeRxnBounds(model,'EX_o2(e)',-1000,'l');
model = changeRxnBounds(model,'ATPM',1000,'u');

% model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% model = addDemandReaction(model,'3pg[c]')
% model = changeObjective(model,'DM_3pg[c]');

% model = addDemandReaction(model,'f6p[c]')
% model = changeObjective(model,'DM_f6p[c]');

model = addDemandReaction(model,'pep[c]')
model = changeObjective(model,'DM_pep[c]');

% model = addDemandReaction(model,'13dpg[c]')
% model = changeObjective(model,'ENO');

FBAsolution = optimizeCbModel(model,'max');

cofactors = {'cmp[c]','ctp[c]','acp[c]','ACP[c]','malACP[c]','amp[c]','atp[c]','adp[c]','pi[c]','ppi[c]','nad[c]','nadh[c]','nadph[c]','nadp[c]','h[c]','h2o[c]','co2[c]'};

% Print fluxes
printFluxVector(model, FBAsolution.x, true)

%Print shadow prices
'Shadow Prices'
printShadowPriceVector(model, FBAsolution.y)

%Print Reduced Costs
'Reduced Costs'
% printFluxVector(model, FBAsolution.w) 