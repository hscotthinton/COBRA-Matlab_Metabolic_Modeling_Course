% FluxVariabilityEthanol.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

% Set objective function
model = changeObjective(model,'EX_etoh(e)');
% model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% % 3 Knockout reactions ( Growth-rate > 0.05)
% model = changeRxnBounds(model,'PTAr',0,'b');
% model = changeRxnBounds(model,'GLUDy',0,'b');
% model = changeRxnBounds(model,'PYK',0,'b');

% Perform flux variability analysis
[minFlux,maxFlux]=fluxVariability(model,100,'max',model.rxns,false,false);

% Print flux values
Difference = abs(maxFlux - minFlux);
FluxDifference = Difference;
n = length(Difference);
for i=1:n
    if Difference(i) < 0.0001
        FluxDifference(i) = 0;
    end
end
% printFluxVector(model, [minFlux, maxFlux, FluxDifference],true)
printFluxVector(model, 1000*[minFlux, maxFlux, FluxDifference]) % print all reactions
