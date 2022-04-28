% AerobicGlucoseBioMass.m

% Script to determine aerobic growth rate of E.coli on glucose
% Taken from "What is flux balance  analysis? - Supplementary tutorial" 
% by J. D. Orth, I. Thiele, & B. 0. Palsson, Nature Biotechnology,
% Volume 28, Number 3, March 2010

clear;

% Input the E.coli core model

model = readCbModel('ecoli_core_model.mat');

% The growth of E. coli on glucose can be simulated under aerobic conditions. 
% To set the maximum glucose uptake rate to 18.5 mmol gDW-1 hr-1 
% (millimoles per gram dry cell weight per hour, the default flux
% units used in the COBRA Toolbox), enter:

model = changeRxnBounds(model,'EX_glc(e)',-18.5,'l');

% This changes the lower bound ('l') of the glucose exchange reaction to 
% -18.5, a biologically realistic uptake rate. By convention, the import 
% of a metabolite is a negative flux. To allow unlimited oxygen uptake, enter:

model = changeRxnBounds(model,'EX_o2(e)',-1000,'l');

% By setting the lower bound of the oxygen uptake reaction to such a 
% large number, it is practically unbounded. 

% Set optimization objective to Biomass_Ecoli_core_N(w/GAM)-Nmet2

model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% Perform FBA with Biomass_Ecoli_core_N(w/GAM)_Nmet2 as the objective, 

FBAsolution = optimizeCbModel(model,'max',0,0) 

% FBAsolution.f then gives the value of the objective function (Z) as 
% 1.6531. This means that the model predicts a growth rate of 1.6531 hr-1. 
% Inspection of the flux distribution vector FBAsolution.x shows that 
% there is high flux in the glycolysis, pentose phosphate, TCA cycle,
% and oxidative phosphorylation pathways, and that no organic by-products 
% are secreted

% Import E.coli core map and adjust parameters

map=readCbMap('ecoli_Textbook_ExportMap');
%options.lb = -10;
%options.ub = 10;
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;
%options.textsize = 12;

% Draw the flux values on the map "target.svg" which can be opened in FireFox

drawFlux(map, model, FBAsolution.x, options);

% Inspection of the flux distribution vector FBAsolution.x shows that 
% there is high flux in the glycolysis, pentose phosphate, TCA cycle,
% and oxidative phosphorylation pathways, and that no organic by-products 
% are secreted. The data values can be seen using the function
%
printFluxVector(model, FBAsolution.x, true)

