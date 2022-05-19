% Figure 3 A, simulate sequential knockdowns of MBNL1 targets in a high
%TGFB + AngII signaling context
%Last update: 3-8-2022 by AN
clc;
clear;
%set figure rendering setting, used to fix a bug encountered when plotting
%default settings may be fine depending on machine
set(0,'DefaultFigureRenderer','Painters')
%declare model file
modelfname='fibroblastMBNL1.xlsx';
% Model indicies for MBNL1 target nodes
sensNode=[35,113,47,28,63,116,119,118,117,114]
%Labels for MBNL1 target nodes
sensRxnNodes={'PDGFR','nYap1',...
'Calcineurin','TGFB1R','p38','SRFmRNA','Runx1','Cbfb','Sox9','MBNL1','TGFB+AngII','Negative Control'};
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
%allocate array for ensemble simulation results
ensemble=[];
COV=0.0331;%input COV for ensemble simulations
%initiate ensemble for loop
for h = [1:150]
ko=0; %ymax value set to simulate knockout

highInput=normrnd(0.6,COV*0.6) %randomly sample increased input values for TGFB and AngII
while highInput < 0 | 1 < highInput %resample if out of range
    highInput=normrnd(0.6,COV*0.6)
end    
%allocate random basal input values
rand_inputs=normrnd(0.1,COV*0.1,[1,9])
while any(rand_inputs < 0 | 1 < rand_inputs) %resample if out of range
    rand_inputs=normrnd(0.1,COV*0.1,[1,9])
end    
%randomly sampled mechanical input
rand_mech=normrnd(0.725,COV*0.725) 
while rand_mech < 0 | 1 < rand_mech  %resample if out of range
    rand_mech=normrnd(0.725,COV*0.725)
end    

%Negative control simulation
% set up simulation params
tspan = [0 500]; % run out to ss
options = [];
%get params
[w,n,EC50,tau,ymax,y0] = paramList{:};
%get random inputs 
w([1,2,4:10])=rand_inputs;
w(3) = rand_mech;
%repack parameters after changing
rpar = [w;n;EC50];
params = {rpar,tau,ymax,specID};
%run model
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yStart = real(y(end,:)'); 

%Positive control simulation
% set up simulation params
tspan = [0 500]; % run out to ss
options = [];
[w,n,EC50,tau,ymax,y0] = paramList{:};
%get random inputs for other model inputs, AngII/TGFB will be overwirrten
%by highInput
w([1,2,4:10])=rand_inputs;
w(3) = rand_mech;
w(2)=highInput; %tgfb input
w(1)=highInput; %angII input
rpar = [w;n;EC50];
params = {rpar,tau,ymax,specID};
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yhighInput = real(y(end,:)'); 



%% Run experimental simulations with MBNL1 target nodes knocked down

ySens=[];

for i = 1:length(sensNode)
% set up simulation params
tspan = [0 500]; % run out to ss
options = [];
%unpack params
[w,n,EC50,tau,ymax,y0] = paramList{:};
%random inputs
w([1,2,4:10])=rand_inputs;
w(3) = rand_mech;
w(2)=highInput; %tgfb input
w(1)=highInput; %angII input
%knock down target node
ymax(sensNode(i))=0;
%repack params
rpar = [w;n;EC50];
params = {rpar,tau,ymax,specID};
%run simulation
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
ySens(:,i) = real(y(end,:)'); 

end
%% get full sensitivity for aSMA
sens=real([ySens,yhighInput,yStart]);
ensemble(h,:)=sens(87,:);
end
%% Save ensemble values
save('KO_Ensembles.mat','ensemble')
%% Load data and Plot Sensitivity
load('KO_Ensembles.mat')
plotNames={'\alpha-SMA'};
fig=figure;
values=mean(ensemble);
errors=(std(ensemble));
[~, I] = sort(values, 'ascend');
b=bar(values(I))
b.FaceColor=[29,145,192]/256
hold on
e1=errorbar(values(I),errors(I))
e1.LineStyle='none';
e1.Color='black'
set(gca, 'XTick', [1:length(sensRxnNodes)])
set(gca, 'XTickLabel', sensRxnNodes(I))
ylabel('Change in \alphaSMA Expression')
xlabel('Knocked-Down Node')
xtickangle(90)

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 10 4]); %x_width=10in y_width=4in
saveas(fig,'Fig3A.png')
saveas(fig,'Fig3A.svg')
