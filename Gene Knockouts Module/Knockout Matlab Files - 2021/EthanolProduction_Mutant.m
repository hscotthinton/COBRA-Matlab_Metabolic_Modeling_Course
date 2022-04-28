% EthanolProduction_Mutant.m

clear;

% Load the E.coli core model
model = readCbModel('ecoli_core_model.mat');

model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

% 5 Knockout reactions
% model = changeRxnBounds(model,'FUM',0,'b');
% model = changeRxnBounds(model,'G6PDH2r',0,'b');
% model = changeRxnBounds(model,'GLUDy',0,'b');
% model = changeRxnBounds(model,'PTAr',0,'b');
% model = changeRxnBounds(model,'SUCDi',0,'b');

% 4 Knockout reactions
% model = changeRxnBounds(model,'ADK1',0,'b');
% model = changeRxnBounds(model,'GLUDy',0,'b');
% model = changeRxnBounds(model,'GND',0,'b');
% model = changeRxnBounds(model,'PTAr',0,'b');

% 3 Knockout reactions ( Growth-rate > 0.14
model = changeRxnBounds(model,'ACKr',0,'b');
model = changeRxnBounds(model,'GLUDy',0,'b');
model = changeRxnBounds(model,'G6PDH2r',0,'b');

% 3 Knockout reactions ( Growth-rate > 0.05
% model = changeRxnBounds(model,'PTAr',0,'b');
% model = changeRxnBounds(model,'GLUDy',0,'b');
% model = changeRxnBounds(model,'PYK',0,'b');

% 2 Knockout reactions
% model = changeRxnBounds(model,'PTAr',0,'b');
% model = changeRxnBounds(model,'G6PDH2r',0,'b');

% 1 Knockout reactions
% model = changeRxnBounds(model,'PFL',0,'b');
% model = changeRxnBounds(model,'ACKr',0,'b');
% model = changeRxnBounds(model,'PTAr',0,'b');

% Set optimization objective
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Perform FBA with Biomass_Ecoli_core_N(w/GAM)_Nmet2 as the objective, 
FBAsolution = optimizeCbModel(model,'max',0,0)
printFluxVector(model,FBAsolution.x, true)

% RobustnessAnalysis
[controlFlux, objFlux] = robustnessAnalysis(model,'EX_etoh(e)',100);

% Import E.coli core map and adjust parameters
map=readCbMap('ecoli_Textbook_ExportMap');
options.lb = -10;
options.ub = 10;
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;

% Draw the flux values on the map "target.svg" which can be opened in FireFox
drawFlux(map, model, FBAsolution.x, options);
