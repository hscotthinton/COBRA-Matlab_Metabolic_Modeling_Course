% Dynamic_SS_Production_BL21_maxProd_100Limit_Plasmid.m
clear;

model=readCbModel('iECD_1391');

% Add plasmid nucleotide precursors
model=addReaction(model,'PLASMID','3276.86 dgtp[c] + 3276.86 dctp[c] + 3199.14 datp[c] + 3199.14 dttp[c] -> plasmid[c]');
model = changeRxnBounds(model,'PLASMID',0.63e-9,'b');

% Add demand reaction for plasmid DNA and set flux rate
model = addDemandReaction(model,'plasmid[c]'); %DM_plasmid[c]

% Add beta-lactamase
model=addReaction(model,'b-lactamase','28 ala-L[c] + 3 cys-L[c] + 16 asp-L[c] + 20 glu-L[c] + 9 phe-L[c] + 21 gly[c] + 7 his-L[c] + 17 ile-L[c] + 11 lys-L[c] + 33 leu-L[c] + 10 met-L[c] + 8 asn-L[c] + 14 pro-L[c] + 9 gln-L[c] + 19 arg-L[c] + 17 ser-L[c] + 20 thr-L[c] + 16 val-L[c] + 4 trp-L[c] + 4 tyr-L[c] + 1,231.52 atp[c] -> b-lactamase[c] + 1,231.52 adp[c] + 1,231.52 pi[c]');
model = changeRxnBounds(model,'b-lactamase',0.000569,'b');

% Add demand reaction for beta-lactamase
model = addDemandReaction(model,'b-lactamase[c]'); % DM_B-lactamase[c]

%Add FlYS3 reaction, change lower bound
model = addReaction(model,'FlYS3','120 ala-L[c] + 4 asp-L[c] + 522 gly[c] + 12 his-L[c] + ile-L[c] + lys-L[c] + 2 met-L[c] + 189 pro-L[c] + 111 ser-L[c] + thr-L[c] + 78 tyr-L[c] + 4476.3 atp[c] + 4476.3 h2o[c] -> flys3[c] + 4476.3 adp[c] + 4476.3 h[c] + 4476.3 pi[c]');

% Add demand reaction
model = addDemandReaction(model,'flys3[c]'); %'DM_flys3[c]'

% Add ferredoxin-NADHP
% model=addReaction(model,'FNADPR','2 flxr[c] + nadp[c] + h[c] <=> nadph[c] + 2 flxso[c]');

% Add ferredoxin-NAD  (EC 1.18.1.3)
% model=addReaction(model,'FNADPR','2 flxr[c] + nad[c] + h[c] <=> nadh[c] + 2 flxso[c]');

model = changeRxnBounds(model, {'EX_glc(e)','EX_o2(e)'},[-10 -20], 'l'); % Aerobic
model = changeRxnBounds(model,'DM_flys3[c]',0,'l'); % Prevent spider silk uptake
model = changeObjective(model,'Ec_biomass_iJO1366_core_53p95M');
model = removeRxns(model,'Ec_biomass_iJO1366_WT_53p95M');

x = 0.0;
for i=1:48
   x = x +.00025;
%    x = 0.0048
   model = changeRxnBounds(model,'FlYS3',x,'b');
   disp(i);
   % Set-up  for dynamicFBA
   substrateRxns = {'EX_glc(e)'};
   initConcentrations = [20000000];
   initBiomass = 47; %  60 g/L wet weight which is approxiamtely 20g/L dry weight
   timeStep = .0625; nSteps = 64; 
   plotRxns = {'EX_ac(e)','DM_flys3[c]','EX_etoh(e)','EX_for(e)','EX_glyclt(e)'};
   % exclUptakeRxns = {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)','EX_nh4(e)','EX_pi(e)','EX_so4(e)','EX_k(e)','EX_fe2(e)',...
   %     'EX_mg2(e)','EX_ca2(e)','EX_cl(e)'};
   exclUptakeRxns = {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)','EX_nh4(e)','EX_pi(e)','EX_so4(e)','EX_k(e)','EX_fe2(e)',...
    'EX_mg2(e)','EX_ca2(e)','EX_cl(e)','EX_cbl1(e)', 'EX_cobalt2(e)', 'EX_cu2(e)', 'EX_fe3(e)' 'EX_meoh(e)', ...
    'EX_mn2(e)', 'EX_mobd(e)', 'EX_na1(e)', 'EX_ni2(e)', 'EX_sel(e)', 'EX_slnt(e)', 'EX_tungs(e)', 'EX_zn2(e)'};

  [concentrationMatrix,excRxnNames,timeVec,biomassVec] = ...
    dynamicFBA(model,substrateRxns,initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns );

  % Plot labels
  subplot(1,2,1); title('Biomass Concentration'); xlabel('Time (h)'); ylabel('Concentration (g/L)');
  subplot(1,2,2); title('Substrate Concentrations'); xlabel('Time (h)'); ylabel('Concentration (mmol/L)');

  % Plot biliverdin production in g/L
  CM = full(concentrationMatrix);
  [a,loc] = ismember('DM_flys3[c]',excRxnNames);
  ssMW = 83.0; % Spider Silk molecular weight (mg/mmol)

  figure(2)
  clf;
  subplot(1,2,1); 
  plot(timeVec,biomassVec);
  title('Biomass Concentration'); xlabel('Time (h)'); ylabel('Concentration (g/L)');
  axis tight;
  subplot(1,2,2); 
  plot(timeVec,ssMW*concentrationMatrix(loc,:));
  % title('Substrate Concentrations'); xlabel('Time (h)'); ylabel('Concentration (mg/L)');
  title('Substrate Concentrations'); xlabel('Time (h)'); ylabel('Concentration (g/L)');
  axis tight;
  legend(strrep(excRxnNames(loc),'EX_',''));

  % Store biliverdin results
  prodFlux(i) = x;
  bioProd(i) = max(ssMW*concentrationMatrix(loc,:));
  ssVec = ssMW*concentrationMatrix(loc,:);
  
  for j=1:nSteps+1
      if biomassVec(j) > 100
          biomassVec(j) = 0;
      end
  end
  [bioMassMax(i), index] = max(biomassVec);
  ssMax(i) = ssVec(index);
  ssTime(i) = index*timeStep;
  testIndex(i) = index;
end

figure(3);
plot(prodFlux, bioProd)
title('Spider Silk Production'); xlabel('Spider Silk (mmol/(gDW h))'); ylabel('Concentration (g/L)');
axis tight

figure(4);
plot(prodFlux, ssMax)
title('Maximum Spider Silk Production at 100g/L Biomass Limit'); xlabel('Spider Silk (mmol/(gDW h))'); ylabel('Spider Silk Concentration (g/L)');
axis tight

figure(5);
plot(prodFlux, ssTime)
title('Production Time to Reach 100 g/L of Biomass'); xlabel('Spider Silk (mmol/(gDW h))'); ylabel('Maximum Production Time (h)');
axis tight

figure(6);
plot(prodFlux, bioMassMax)
title('Maximum Biomass Production'); xlabel('Spider Silk (mmol/(gDW h))'); ylabel('Biomass Concentration (g/L)');
axis tight