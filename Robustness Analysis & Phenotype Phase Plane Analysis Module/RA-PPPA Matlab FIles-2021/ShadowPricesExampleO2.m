% ShadowPricesExampleO2.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

model = changeRxnBounds(model,'EX_glc(e)',-10,'b');

growthRates = zeros(61,1);
for i = 0:60
    model = changeRxnBounds(model,'EX_o2(e)',-i,'b');
    FBAsolution = optimizeCbModel(model,'max');
    growthRates(i+1) = FBAsolution.f;
    shadowPrice(i+1) = - FBAsolution.y(57); % o2[e] is metabolite number 57; add negaitve sign for gurobi solver
    xaxis1(i+1) = -i;
end

figure(1)
[haxes,hline1,hline2]=plotyy(xaxis1,growthRates,xaxis1,shadowPrice);
ylabel(haxes(1),'BioMass Objective Function')
ylabel(haxes(2),'EX-o2(e) Shadow Prices')
xlabel('EX-o2(e)');