
%% Runs the function 'validateModel_MBNL1_Ensemble'
% Last updated by AN on 03-10-2022
clc;
clear;
%name of validation file
validationfname='fib_validationInt_MBNL1';
%model needs to be in xlsx format 
modelfname='fibroblastMBNL1'; 
%run validation function
COV=0.0331; %input COV
[percentMatch, resultChart] = validateModel_MBNL1_Ensemble(modelfname, validationfname, COV)
%% Clear the old results chart and rewrite
if exist('MBNL1_resultsChart_Int.csv','file')==2
  delete('MBNL1_resultsChart_Int.csv');
end

xlswrite('MBNL1_resultsChart_Int.csv',resultChart)
save('MBNL1_intResult.mat','resultChart')
%% Plot figure 2C
load('MBNL1_intResult.mat')
usubNetwork=unique(subNetwork);
totals=[];
matched=[];
for k = 1:length(usubNetwork)
    tCount=0; %# total relationships
    mCount=0; %# relationaships matching
    for j = 1:length(match)
        if strcmp(subNetwork{j},usubNetwork{k})
            tCount=tCount+1;
            if strcmp(match{j},'yes')
                mCount=mCount+1;
            end
        end

    end
    totals(k)=tCount;
    matched(k)=mCount;
end
        
figure;
x=([1:length(usubNetwork)]);
values=[matched'./totals']*100;
b=bar(values)
b.FaceColor=[29,145,192]/256
