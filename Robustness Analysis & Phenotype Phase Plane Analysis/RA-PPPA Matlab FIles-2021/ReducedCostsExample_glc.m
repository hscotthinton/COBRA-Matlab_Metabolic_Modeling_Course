% ReducedCostsExample_glc.m
clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');
model = changeRxnBounds(model,'EX_o2(e)',-20,'b');

growthRates = zeros(30,1);

growthRates(1) = 0;
reducedCosts(1) = 0;
xaxis1(1) = -1;

growthRates(2) = 0;
reducedCosts(2) = 0;
xaxis1(2) = -2;

growthRates(3) = 0;
reducedCosts(3) = 0;
xaxis1(3) = -3;

for i = 4:30
    model = changeRxnBounds(model,'EX_glc(e)',-i,'b');
    FBAsolution = optimizeCbModel(model,'max');
    growthRates(i) = FBAsolution.f;
    reducedCosts(i) = - FBAsolution.w(28); % EX_glc(e) is reaction number 35
    xaxis1(i) = -i;
end

figure(1)
[haxes,hline1,hline2]=plotyy(xaxis1,growthRates,xaxis1,reducedCosts);
ylabel(haxes(1),'BioMass Objective Function')
ylabel(haxes(2),'EX-glc(e) Reduced Costs')
xlabel('EX-glc(e)');