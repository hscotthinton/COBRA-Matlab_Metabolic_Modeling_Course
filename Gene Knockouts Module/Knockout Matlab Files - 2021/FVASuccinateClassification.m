% FVASuccinateClassification.m
clear; clc;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Change carbon source from glucose to succinate
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)_Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Perform flux variability analysis classification
[minFlux,maxFlux]=fluxVariability(model,90,'max',model.rxns,false,false);

% Print flux values
Difference = abs(maxFlux - minFlux);
FluxDifference = Difference;
n = length(Difference);
for i=1:n
    if Difference(i) < 0.0001
        FluxDifference(i) = 0;
    end
end
averageFlux = (maxFlux - minFlux)/2;
printFluxVector(model, [minFlux, maxFlux, averageFlux, FluxDifference])

BioMassID = findRxnIDs(model, 'Biomass_Ecoli_core_N(w/GAM)-Nmet2');
BiomassRatio = minFlux(BioMassID)/maxFlux(BioMassID);

% Find hard-coupled reactions
HCReactions = {};
n = length(maxFlux);
j = 1;
for i=1:n
    if (minFlux(i) == BiomassRatio*maxFlux(i)) && (maxFlux(i) > 0)
        HCReactions(j) = model.rxns(i);
        j = j+1;
    end
end
HardCoupledReactions = transpose(HCReactions)

% Find partially-coupled reactions
PCReactions = {};
n = length(maxFlux);
j = 1;
for i=1:n
    if (minFlux(i) > 0 ) && (minFlux(i) < BiomassRatio*maxFlux(i))
        PCReactions(j) = model.rxns(i);
        j = j+1;
    end
end
PartiallyCoupledReactions = transpose(PCReactions)

% Find not-coupled reactions
NCReactions = {};
n = length(maxFlux);
j = 1;
for i=1:n
    if (minFlux(i) <= 0 ) && (minFlux(i) < maxFlux(i))
        NCReactions(j) = model.rxns(i);
        j = j+1;
    end
end
NotCoupledReactions = transpose(NCReactions)

% Find zero-flux reactions
ZFReactions = {};
n = length(maxFlux);
j = 1;
for i=1:n
    if (minFlux(i) == 0 ) && (maxFlux(i) == 0)
        ZFReactions(j) = model.rxns(i);
        j = j+1;
    end
end
ZeroFluxReactions = transpose(ZFReactions)

