
%% Runs the function 'validateModel_IO_Ensemble'
%Runs Input-Out validation function
%Last updated by AN on 03-10-2022
clc;
clear all;
%name of validation file
validationfname='fib_validationIO_MBNL1';
%name of model file
modelfname='fibroblastMBNL1'; 
controlInputW=0.1; %Basal input mean
COV=0.0331; %input COV for random sampling
%run validation function
[percentMatch, resultChart,MatchControl] =...
     validateModel_IO_Ensemble(modelfname, validationfname,controlInputW,COV);
%% Print full results of the validation
if exist('MBNL1_resultsChart_IO.csv','file')==2
  delete('MBNL1_resultsChart_IO.csv');
end

xlswrite('MBNL1_resultsChart_IO.csv',resultChart)

