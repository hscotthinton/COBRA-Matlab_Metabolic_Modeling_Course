% GlucoseEthanolOptknock_iaf1260.m
clear; 

% Set constraints
model = readCbModel('iAF1260.mat');
model = changeRxnBounds(model, {'EX_o2(e)', 'EX_glc(e)'}, [0 -10], 'l');
model = changeObjective(model,'Ec_biomass_iAF1260_core_59p81M');

% Remove exchange & transport reactions & ATPM & Biomass reaction
[transRxns,nonTransRxns] = findTransRxns(model,true); % Remove transport reactions
[tmp,ATPMnumber] = ismember('ATPM',nonTransRxns); % Identify ATPM reaction number
nonTransRxns(ATPMnumber) = []; % Remove ATPM from potential knockout reactions
selectedRxns = nonTransRxns;

% Complete pFBA to remove essential and zero flux reactions
disp('Executing pFBA');
[GeneClasses RxnClasses modelIrrevFM] = pFBA(model, 'geneoption',0, 'tol',1e-7);
unwantedReactions = [RxnClasses.Essential_Rxns;RxnClasses.ZeroFlux_Rxns;RxnClasses.Blocked_Rxns];
[tmp,BioMassnumber] = ismember('Ec_biomass_iAF1260_core_59p81M',unwantedReactions); % Identify biomass reaction number
unwantedReactions(BioMassnumber) = []; % Remove Biomass from potential knockout reactions

% remove essential and zero flux reactions from selected reactions pool
[tmp,IRNumbers] = ismember(unwantedReactions,selectedRxns); 

j = 1;
for i=1:length(IRNumbers)
    if IRNumbers(i) > 0
        NZNumbers(j) = IRNumbers(i);
        j = j + 1;
    end
end
selectedRxns(NZNumbers') = [];

% Optknock analysis for Ethanol secretion
disp('Executing optKnock');
options.targetRxn = 'EX_etoh(e)';
options.vMax = 1000;
options.numDel = 3;
options.numDelSense = 'L';
constrOpt.rxnList = {'Ec_biomass_iAF1260_core_59p81M','ATPM'};
constrOpt.values = [0.05, 8.39];
constrOpt.sense = 'GE';

optKnockSol = OptKnock(model, selectedRxns, options, constrOpt);
deletions = optKnockSol.rxnList'

% Print out growth rate and minimum & maximum secretion rate
[growthRate,minProd,maxProd] = testOptKnockSol(model,'EX_etoh(e)',optKnockSol.rxnList)

% Print production envelope
lineColor = 'b';
targetRxn = 'EX_etoh(e)';
biomassRxn = 'Ec_biomass_iAF1260_core_59p81M';
geneDelFlag = false;
nPts = 50;
figure(1)
[biomassValues,targetValues] = productionEnvelope(model,deletions,lineColor,targetRxn,biomassRxn,geneDelFlag,nPts);
xlabel('Biomass (mmol/g DW-hr)')
ylabel('Ethanol(mmol/g DW-hr)')

figure(2)
[biomassValues,targetValues] = multiProductionEnvelope(model,deletions,biomassRxn,false,50,false);

