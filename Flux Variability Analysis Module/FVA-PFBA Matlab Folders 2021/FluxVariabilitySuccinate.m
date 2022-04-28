% FluxVariabilitySuccinate.m
clear; clc;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Change carbon source from glucose to succinate
model = changeRxnBounds(model,'EX_glc(e)',0,'l');
model = changeRxnBounds(model,'EX_succ(e)',-20,'l');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Perform flux variability analysis

% [minFlux,maxFlux]=fluxVariability(model);
[minFlux,maxFlux]=fluxVariability(model,100,'max',model.rxns,false,false);
% If loops are not allowed then the minFlux and maxFlux values will need to
% be multiplied by 1000

% Print flux values
Difference = abs(maxFlux - minFlux);
FluxDifference = Difference;
n = length(Difference);
for i=1:n
    if Difference(i) < 0.0001
        FluxDifference(i) = 0;
    end
end
printFluxVector(model, [minFlux, maxFlux, FluxDifference], true)

% Save flux variability map

map=readCbMap('ecoli_Textbook_ExportMap');
drawFluxVariability(map,model,minFlux,maxFlux);