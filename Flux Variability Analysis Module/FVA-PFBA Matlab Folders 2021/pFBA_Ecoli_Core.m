% pFBA_Ecoli_Core2.m
clear;

model = readCbModel('ecoli_core_model.mat');

model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
model = changeRxnBounds(model,'EX_o2(e)',-30,'l');

map=readCbMap('ecoli_Textbook_ExportMap');

% [GeneClasses RxnClasses modelIrrevFM] = pFBA(model, 'geneoption',0,'map','ecoli_Textbook_ExportMap.txt','tol',1e-7)
[GeneClasses RxnClasses modelIrrevFM] = pFBA(model, 'geneoption',0,'tol',1e-7) % No map