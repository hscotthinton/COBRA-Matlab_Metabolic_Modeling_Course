% AerobicGlucoseBioMassPPP.m

clear;

% Input the E.coli core model
model = readCbModel('ecoli_core_model.mat');

% Set the lower bounds for oxygen and glucose uptake
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
model = changeRxnBounds(model,'EX_glc(e)',-20,'l');

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Using phenotype plane analysis, plot the objective function as a 
% function of the glucose and oxygen uptake rate
[growthRates,shadowPrices1,shadowPrices2]=phenotypePhasePlane(model,'EX_glc(e)','EX_o2(e)');

% Prep for figures
nPts = 50;
range1 = 20;
range2 = 20;
ind1 = linspace(0,range1,nPts);
ind2 = linspace(0,range2,nPts);

% Plot shadowPrices1 (Glucose) isoclines
figure(4);
surf(ind1,ind2,-shadowPrices1);
xlabel('EX glc(e) (mmol/g DW-hr)'); 
ylabel('EX o2(e) (mmol/g DW-hr)'); 
% axis([0 20 0 20 0 0.15]);
title('Glucose Shadow Prices (shadowPrices1)');

% Plot shadowPrices2 (Oxygen) isoclines
figure(5);
surf(ind1,ind2,-shadowPrices2);
xlabel('EX glc(e) (mmol/g DW-hr)'); 
ylabel('EX o2(e) (mmol/g DW-hr)'); 
title('Oxygen Shadow Prices (shadowPrices2)');

% Plot alpha isoclines
figure(6);
alpha = - shadowPrices1./shadowPrices2;
pcolor(ind1,ind2,alpha);
xlabel('EX glc(e) (mmol/g DW-hr)'); 
ylabel('EX o2(e) (mmol/g DW-hr)'); 
title('Ratio of Shadow Prices (alpha)');
colorbar;

figure(7);
surf(ind1,ind2,alpha);
xlabel('EX glc(e) (mmol/g DW-hr)'); 
ylabel('EX o2(e) (mmol/g DW-hr)'); 
zlabel('Ratio of Shadow Prices'); 
title('Ratio of Shadow Prices (alpha)');

figure(8);
contour(ind1,ind2,growthRates,20);
xlabel('EX glc(e) (mmol/g DW-hr)'); 
ylabel('EX o2(e) (mmol/g DW-hr)'); 
zlabel('Growth-rate'); 
title('Phenotype Phase Plane Contour Plot');

figure(9);
pcolor(ind1,ind2,-shadowPrices1);
xlabel('EX glc(e) (mmol/g DW-hr)'); 
ylabel('EX o2(e) (mmol/g DW-hr)'); 
zlabel('Growth-rate'); 
title('Glucose Shadow Prices (shadowPrices1)');

figure(10);
pcolor(ind1,ind2,-shadowPrices2);
xlabel('EX glc(e) (mmol/g DW-hr)'); 
ylabel('EX o2(e) (mmol/g DW-hr)'); 
zlabel('Growth-rate'); 
title('Oxygen Shadow Prices (shadowPrices2)');
