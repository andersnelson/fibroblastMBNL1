
%% Runs the function 'ensembleConvergance'
%Plots results of the function
%Last updated by AN on 03-10-2022
clc;
clear;
%name of validation file
validationfname='fib_validationInt_MBNL1.xlsx'
%model needs to be in xlsx format 
modelfname='fibroblastMBNL1.xlsx' 
%run validation function
COV=0.0331;%input COV
[meanDelta5,validationDeltas] = ensembleConvergence(modelfname, validationfname,COV); %call convergence function

%% Save simulation results from convergance screen
save('converganceEnsemble.mat','meanDelta5')
%% Find index for best N
N = find(meanDelta5<0.005)
%% plot results
 load('converganceEnsemble.mat')

fig=figure;
plot([3:200],validationDeltas(3:200))
title('Ensemble Simualtion Convergence Individual %Matching')
xlabel('Number of Simulations')
ylabel('%Validation')
hold on
plot(150,validationDeltas(150),'b*')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 10 4]); %x_width=9in y_width=4in
saveas(fig,'FigS2.svg')
saveas(fig,'FigS2.png')
