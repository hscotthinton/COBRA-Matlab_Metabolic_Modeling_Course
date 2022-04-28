% ReducedCostsExampleO2.m

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
    reducedCosts(i+1) = - FBAsolution.w(36); % EX_o2(e) is metabolite number 36; Change sign for consistency
    xaxis1(i+1) = -i;
end

figure(1)
[haxes,hline1,hline2]=plotyy(xaxis1,growthRates,xaxis1,reducedCosts);
xlabel('EX-o2(e)');
axes(haxes(1));
axis([-60 0 0 1]);
ylabel(haxes(1),'BioMass Objective Function');
axes(haxes(2));
axis([-60 0 -0.06 0.03]);
ylabel(haxes(2),'EX-o2(e) Reduced Costs');
xlabel('EX-o2(e)');