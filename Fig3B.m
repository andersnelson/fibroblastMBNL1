% Figure 3 B, simulate sequential overexpressions of MBNL1 targets and
% MBNL1
%Last update: 3-8-2022 by AN
clc;
clear;
%set figure rendering setting, used to fix a bug encountered when plotting
%default settings may be fine depending on machine
set(0,'DefaultFigureRenderer','Painters')
%declare model file
modelfname='fibroblastMBNL1.xlsx';
%model indicies for MBNL1 target nodes
sensNode=[35,113,47,28,63,116,119,118,117,114];
%Labels for MBNL1 target nodes
sensRxnNodes={'PDGFR','nYap1','Calcineurin','TGFB1R','p38',... %labels for target nodes, used later for plotting
    'SRFmRNA','Runx1','Cbfb','Sox9','MBNL1','Negative Control'};
% delete the previously formed ODE if it exists and generate new temp model 
% ODE file
pwd = cd;
if exist([pwd '\ODEfun.m'],'file') == 2
delete('ODEfun.m');
end
% parse out model name (xls2Netflux needs it as an arg)
namepos = findstr('.xls', modelfname);
namestr = modelfname(1:namepos-1);
namestr = cellstr(namestr);
% generate ODE from model spreadsheet
[specID,reactionIDs,~,paramList,ODElist,~, error] = util.xls2Netflux(namestr,modelfname);
commandLine = util.exportODE2(specID,paramList,ODElist);
util.textwrite('ODEfun.m',commandLine);
%set random seed
randn('seed',0);
%allocate array for ensemble simulations
ensemble=[];
COV=0.0331; %Input COV
%%
for k = [1:150]
%allocate random basal input values
rand_inputs=normrnd(0.1,COV*0.1,[1,9])
while any(rand_inputs < 0 | 1 < rand_inputs) %resample if out of range
    rand_inputs=normrnd(0.1,COV*0.1,[1,9])
end
%randomly sample mech. stimulus
rand_mech=normrnd(0.725,0.725*COV)
while rand_mech < 0 | 1 < rand_mech 
    rand_mech=normrnd(0.725,0.725*COV) %resample if out of range
end    
% set up simulation params
tspan = [0 500]; % run out to ss
options = [];
[w,n,EC50,tau,ymax,y0] = paramList{:};
%get baseline simulation results with no overexpression
w([1,2,4:10])=rand_inputs;
w(3)=rand_mech;
rpar = [w;n;EC50];
params = {rpar,tau,ymax,specID};
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yStart = real(y(end,:)'); 


%% Run simulations with MBNL1 target nodes OE

ySens=[];
for i = 1:length(sensNode)
% set up simulation params
tspan = [0 500]; % run out to ss
options = [];
%unpack parameters and alter
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2,4:10])=rand_inputs;
w(3)=rand_mech;
over=1.0; %overexpression
y0(sensNode(i))=over;
ymax(sensNode(i))=over;
tau(sensNode(i))=100000;
%repack parameters for simulation
rpar = [w;n;EC50];
params = {rpar,tau,ymax,specID};
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
ySens(:,i) = real(y(end,:)'); 

end
%% get full sensitivity for aSMA
sens=real([ySens,yStart]);
ensemble(k,:)=sens(87,:);
end
%% Save ensemble values
save('OE_Ensembles_Fig3.mat','ensemble')
%% Load Ensemble values and Plot Sensitivity
load('OE_Ensembles_Fig3.mat')
plotNames={'\alpha-SMA'};
fig=figure;
values=mean(ensemble)
errors=std(ensemble);
[~, I] = sort(values, 'ascend');
b=bar(values(I))
b.FaceColor=[29,145,192]/256
hold on
e1=errorbar(values(I),errors(I))
e1.LineStyle='none';
e1.Color='black'
set(gca, 'XTick', [1:length(sensRxnNodes)])
set(gca, 'XTickLabel', sensRxnNodes(I))
ylabel('\alphaSMA Expression')
xlabel('Overexpressed Node')
xtickangle(90)
ylim([0,1.1])
title(['Overexpression Screen'])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 10 4]); %x_width=10in y_width=4in
saveas(fig,'Fig3B.png')
saveas(fig,'Fig3B.svg')


