% Dynamic_Growth_K12media_iJO1366_FBC.m
clear;

% Load iJO1366 model
model = readCbModel('iJO1366.mat');
model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');

%Setting carbon source and oxygen
model = changeRxnBounds(model,'EX_glc__D_e',-10,'l');
model = changeRxnBounds(model,'EX_o2_e',-30,'l');

% Set uptake values for amino acids;
model = changeRxnBounds(model,'EX_ala__L_e',-0.200200247,'l');
model = changeRxnBounds(model,'EX_arg__L_e',-0.065982824,'l');
model = changeRxnBounds(model,'EX_asn__L_e',-0.095998062,'l');
model = changeRxnBounds(model,'EX_asp__L_e',-0.09529124,'l');
model = changeRxnBounds(model,'EX_cys__L_e',-0.013085243,'l');
model = changeRxnBounds(model,'EX_gln__L_e',-0.187137594,'l');
model = changeRxnBounds(model,'EX_glu__L_e',-0.185878393,'l');
model = changeRxnBounds(model,'EX_gly_e',-0.14255367,'l');
model = changeRxnBounds(model,'EX_his__L_e',-0.030655649,'l');
model = changeRxnBounds(model,'EX_ile__L_e',-0.08762833,'l');
model = changeRxnBounds(model,'EX_leu__L_e',-0.126909995,'l');
model = changeRxnBounds(model,'EX_lys__L_e',-0.119293303,'l');
model = changeRxnBounds(model,'EX_met__L_e',-0.02390703,'l');
model = changeRxnBounds(model,'EX_phe__L_e',-0.05758489,'l');
model = changeRxnBounds(model,'EX_pro__L_e',-0.079180891,'l');
model = changeRxnBounds(model,'EX_ser__L_e',-0.098060253,'l');
model = changeRxnBounds(model,'EX_thr__L_e',-0.089838012,'l');
model = changeRxnBounds(model,'EX_tyr__L_e',-0.039374888,'l');
model = changeRxnBounds(model,'EX_val__L_e',-0.111648451,'l');

% Set uptake values for minerals;
model = changeRxnBounds(model,'EX_ca2_e',-0.00237348,'l');
model = changeRxnBounds(model,'EX_cobalt2_e',-1.14e-05,'l');
model = changeRxnBounds(model,'EX_ni2_e',-0.000147288,'l');
model = changeRxnBounds(model,'EX_nh4_e',-6.002150375,'l');
model = changeRxnBounds(model,'EX_pi_e',-5.555386948,'l');
model = changeRxnBounds(model,'EX_k_e',-3.943619038,'l');
model = changeRxnBounds(model,'EX_so4_e',-0.167768684,'l');
model = changeRxnBounds(model,'EX_cl_e',-0.138740225,'l');
model = changeRxnBounds(model,'EX_cu2_e',-0.000634963,'l');
model = changeRxnBounds(model,'EX_fe3_e',-0.00696512,'l');
model = changeRxnBounds(model,'EX_mg2_e',-0.160804933,'l');
model = changeRxnBounds(model,'EX_mn2_e',-0.008010947,'l');
model = changeRxnBounds(model,'EX_mobd_e',-0.001508618,'l');
model = changeRxnBounds(model,'EX_na1_e',-0.102668246,'l');
model = changeRxnBounds(model,'EX_thm_e',-0.000746848,'l');
model = changeRxnBounds(model,'EX_zn2_e',-0.001378378,'l');
 
% Set-up variables for dynamicFBA
% NOTE- substrate rxns and plot rxns need to be in the order that they
% appear in the model
 
initBiomass = .01;
timeStep = 0.5; nSteps = 100;
substrateRxns = {'EX_ala__L_e','EX_arg__L_e','EX_asn__L_e','EX_asp__L_e','EX_cl_e','EX_cu2_e','EX_cys__L_e','EX_fe3_e','EX_glc__D_e','EX_gln__L_e','EX_glu__L_e','EX_gly_e','EX_his__L_e','EX_ile__L_e','EX_k_e','EX_leu__L_e','EX_lys__L_e','EX_met__L_e','EX_mg2_e','EX_mn2_e','EX_mobd_e','EX_na1_e','EX_nh4_e','EX_phe__L_e','EX_pi_e','EX_pro__L_e','EX_ser__L_e','EX_so4_e','EX_thm_e','EX_thr__L_e','EX_tyr__L_e','EX_val__L_e','EX_zn2_e'};
initConcentrations = [2.525535975,0.832376579,1.211020285,1.202103681,1.750214773,0.008010093,0.165070981,0.087865335,138.7655417,2.360749966,2.344865085,1.798321567,0.386722527,1.105435694,49.74894838,1.600975833,1.504890895,0.301588365,2.028562155,0.101058487,0.019031286,1.295164976,75.71742258,0.726436225,70.08147994,0.998870842,1.237034922,2.11641021,0.009421519,1.133310947,0.496716154,1.408450704,0.017388304];
plotRxns = {'EX_ac_e','EX_acald_e','EX_ala__L_e','EX_arg__L_e','EX_asn__L_e','EX_asp__L_e','EX_cl_e','EX_cu2_e','EX_cys__L_e','EX_etoh_e','EX_fe3_e','EX_for_e','EX_glc__D_e','EX_gln__L_e','EX_glu__L_e','EX_gly_e','EX_his__L_e','EX_ile__L_e','EX_k_e','EX_lac__L_e','EX_leu__L_e','EX_lys__L_e','EX_met__L_e','EX_mg2_e','EX_mn2_e','EX_mobd_e','EX_na1_e','EX_nh4_e','EX_phe__L_e','EX_pi_e','EX_pro__L_e','EX_ser__L_e','EX_so4_e','EX_succ_e','EX_thm_e','EX_thr__L_e','EX_tyr__L_e','EX_val__L_e','EX_zn2_e'};
 
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
