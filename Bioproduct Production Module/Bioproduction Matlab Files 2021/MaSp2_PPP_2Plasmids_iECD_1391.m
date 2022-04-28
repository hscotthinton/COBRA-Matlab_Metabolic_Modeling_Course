% MaSp2_PPP_2Plasmids_iECD_1391.m
clear; 

% Input the E.coli core model
model=readCbModel('iECD_1391.xml');

% Add plasmid nucleotide precursors
model=addReaction(model,'PLASMID','3275.84 dgtp[c] + 3275.84 dctp[c] + 3198.12 datp[c] + 3198.12 dttp[c] -> plasmid[c]');
model = changeRxnBounds(model,'PLASMID',0.63e-9,'l');
model = addDemandReaction(model,'plasmid[c]');

% Add kanamycin resistance 
model=addReaction(model,'KanR','35 ala-L[c] + 19 arg-L[c] + 3 asn-L[c] + 25 asp-L[c] + 5 cys-L[c] + 11 gln-L[c] + 18 glu-L[c] + 21 gly[c] + 7 his-L[c] + 10 ile-L[c] + 32 leu-L[c] + 4 lys-L[c] + 6 met-L[c] + 11 phe-L[c] + 11 pro-L[c] + 11 ser-L[c] + 10 thr-L[c] + 5 trp-L[c] + 4 tyr-L[c] + 16 val-L[c] + 1136.8 atp[c] + 1136.8 h2o[c] -> kanr[c] + 1136.8 adp[c] + 1136.8 h[c] + 1136.8 pi[c]');
model = changeRxnBounds(model,'KanR',0.000569,'l');
model = addDemandReaction(model,'kanr[c]');

% Add plasmid2 nucleotide precursors reactions (4GPP)
model = addReaction(model,'PLASMID2','2593.76 dgtp[c] + 2593.76 dctp[c] + 2532.24 datp[c] + 2532.24 dttp[c] -> plasmid2[c]');
model = changeRxnBounds(model,'PLASMID2',0.63e-9,'l');
model = addDemandReaction(model,'plasmid2[c]');

% Add Chloramphenicol resistance 
model = addReaction(model,'CmR','15 ala-L[c] + 5 arg-L[c] + 10 asn-L[c] + 12 asp-L[c] + 5 cys-L[c] + 13 gln-L[c] + 12 glu-L[c] + 10 gly[c] + 12 his-L[c] + 9 ile-L[c] + 13 leu-L[c] + 12 lys-L[c] + 9 met-L[c] + 20 phe-L[c] + 7 pro-L[c] + 10 ser-L[c] + 13 thr-L[c] + 5 trp-L[c] + 11 tyr-L[c] + 16 val-L[c] + 1136.8 atp[c] + 1136.8 h2o[c] -> cmr[c] + 943 adp[c] + 943 h[c] + 943 pi[c]');
model = changeRxnBounds(model,'CmR',0.000643,'l');
model = addDemandReaction(model,'cmr[c]');

% Add MaSp2 reaction
model = addReaction(model,'MaSp2','288 ala-L[c] + 4 asp-L[c] + 192 gln-L[c] + 512 gly[c] + 4 his-L[c] + 4 lys-L[c] + 4 met-L[c] + 236 pro-L[c] + 76 ser-L[c] + 4 thr-L[c] + 68 tyr-L[c] + 6011.12 atp[c] + 6011.12 h2o[c] -> masp2[c] + 6011.12 adp[c] + 6011.12 h[c] + 6011.12 pi[c]');

% Add demand reaction
model = addDemandReaction(model,'masp2[c]'); %'DM_masp2[c]'
% model = changeRxnBounds(model,'MaSp2',0.00000,'l'); % MaSp2 < 0.00233, 0.006, 0.008 for EX_glc(e) > -10
model = changeRxnBounds(model,'MaSp2',0.00233,'l'); % MaSp2 < 0.001, 0.00233, 0.0038 for EX_glc(e) > -5

% Set key variables
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
model = changeObjective(model,'Ec_biomass_iJO1366_core_53p95M');

% Phenotype Phase Plane Analysis
[growthRates,shadowPrices1,shadowPrices2]=phenotypePhasePlane(model,'EX_glc(e)','EX_o2(e)');

title('Growth-rate Phenotype Phase Plane Analysis with Two Plasmids'); 
zlabel('Growth-rate(1/hr)');

% title('Growth Rate Phenotype Phase Plane Analysis with No Plasmids'); 
