% Dynamic_Growth_K12media_iJO1366.m
clear;

% Load iJO1366 model
load('ecoli_iJO1366');
model = changeObjective(model,'Ec_biomass_iJO1366_core_53p95M');

%Setting carbon source and oxygen
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-0,'l');

% Set uptake values for amino acids;
model = changeRxnBounds(model,'EX_ala-L(e)',-0.200200247,'l');
model = changeRxnBounds(model,'EX_arg-L(e)',-0.065982824,'l');
model = changeRxnBounds(model,'EX_asn-L(e)',-0.095998062,'l');
model = changeRxnBounds(model,'EX_asp-L(e)',-0.09529124,'l');
model = changeRxnBounds(model,'EX_cys-L(e)',-0.013085243,'l');
model = changeRxnBounds(model,'EX_gln-L(e)',-0.187137594,'l');
model = changeRxnBounds(model,'EX_glu-L(e)',-0.185878393,'l');
model = changeRxnBounds(model,'EX_gly(e)',-0.14255367,'l');
model = changeRxnBounds(model,'EX_his-L(e)',-0.030655649,'l');
model = changeRxnBounds(model,'EX_ile-L(e)',-0.08762833,'l');
model = changeRxnBounds(model,'EX_leu-L(e)',-0.126909995,'l');
model = changeRxnBounds(model,'EX_lys-L(e)',-0.119293303,'l');
model = changeRxnBounds(model,'EX_met-L(e)',-0.02390703,'l');
model = changeRxnBounds(model,'EX_phe-L(e)',-0.05758489,'l');
model = changeRxnBounds(model,'EX_pro-L(e)',-0.079180891,'l');
model = changeRxnBounds(model,'EX_ser-L(e)',-0.098060253,'l');
model = changeRxnBounds(model,'EX_thr-L(e)',-0.089838012,'l');
model = changeRxnBounds(model,'EX_tyr-L(e)',-0.039374888,'l');
model = changeRxnBounds(model,'EX_val-L(e)',-0.111648451,'l');

% Set uptake values for minerals;
model = changeRxnBounds(model,'EX_ca2(e)',-0.00237348,'l');
model = changeRxnBounds(model,'EX_cobalt2(e)',-1.14e-05,'l');
model = changeRxnBounds(model,'EX_ni2(e)',-0.000147288,'l');
model = changeRxnBounds(model,'EX_nh4(e)',-6.002150375,'l');
model = changeRxnBounds(model,'EX_pi(e)',-5.555386948,'l');
model = changeRxnBounds(model,'EX_k(e)',-3.943619038,'l');
model = changeRxnBounds(model,'EX_so4(e)',-0.167768684,'l');
model = changeRxnBounds(model,'EX_cl(e)',-0.138740225,'l');
model = changeRxnBounds(model,'EX_cu2(e)',-0.000634963,'l');
model = changeRxnBounds(model,'EX_fe3(e)',-0.00696512,'l');
model = changeRxnBounds(model,'EX_mg2(e)',-0.160804933,'l');
model = changeRxnBounds(model,'EX_mn2(e)',-0.008010947,'l');
model = changeRxnBounds(model,'EX_mobd(e)',-0.001508618,'l');
model = changeRxnBounds(model,'EX_na1(e)',-0.102668246,'l');
model = changeRxnBounds(model,'EX_thm(e)',-0.000746848,'l');
model = changeRxnBounds(model,'EX_zn2(e)',-0.001378378,'l');
 
% Set-up variables for dynamicFBA
% NOTE- substrate rxns and plot rxns need to be in the order that they
% appear in the model
 
initBiomass = .01;
timeStep = 0.5; nSteps = 100;
substrateRxns = {'EX_ala-L(e)','EX_arg-L(e)','EX_asn-L(e)','EX_asp-L(e)','EX_cl(e)','EX_cu2(e)','EX_cys-L(e)','EX_fe3(e)','EX_glc(e)','EX_gln-L(e)','EX_glu-L(e)','EX_gly(e)','EX_his-L(e)','EX_ile-L(e)','EX_k(e)','EX_leu-L(e)','EX_lys-L(e)','EX_met-L(e)','EX_mg2(e)','EX_mn2(e)','EX_mobd(e)','EX_na1(e)','EX_nh4(e)','EX_phe-L(e)','EX_pi(e)','EX_pro-L(e)','EX_ser-L(e)','EX_so4(e)','EX_thm(e)','EX_thr-L(e)','EX_tyr-L(e)','EX_val-L(e)','EX_zn2(e)'};
initConcentrations = [2.525535975,0.832376579,1.211020285,1.202103681,1.750214773,0.008010093,0.165070981,0.087865335,138.7655417,2.360749966,2.344865085,1.798321567,0.386722527,1.105435694,49.74894838,1.600975833,1.504890895,0.301588365,2.028562155,0.101058487,0.019031286,1.295164976,75.71742258,0.726436225,70.08147994,0.998870842,1.237034922,2.11641021,0.009421519,1.133310947,0.496716154,1.408450704,0.017388304];
plotRxns = {'EX_ac(e)','EX_acald(e)','EX_ala-L(e)','EX_arg-L(e)','EX_asn-L(e)','EX_asp-L(e)','EX_cl(e)','EX_cu2(e)','EX_cys-L(e)','EX_etoh(e)','EX_fe3(e)','EX_for(e)','EX_glc(e)','EX_gln-L(e)','EX_glu-L(e)','EX_gly(e)','EX_his-L(e)','EX_ile-L(e)','EX_k(e)','EX_lac-L(e)','EX_leu-L(e)','EX_lys-L(e)','EX_met-L(e)','EX_mg2(e)','EX_mn2(e)','EX_mobd(e)','EX_na1(e)','EX_nh4(e)','EX_phe-L(e)','EX_pi(e)','EX_pro-L(e)','EX_ser-L(e)','EX_so4(e)','EX_succ(e)','EX_thm(e)','EX_thr-L(e)','EX_tyr-L(e)','EX_val-L(e)','EX_zn2(e)'};
 
dynamicFBA(model,substrateRxns,initConcentrations, initBiomass, timeStep, nSteps, plotRxns);

%labeling
subplot(1,2,1);
title('Biomass Concentration');
xlabel('Time steps');
ylabel('Concentration (g/L)');
subplot(1,2,2);
title('Substrate Concentrations');
xlabel('Time steps');
ylabel('Concentrations (mmol/L)');
