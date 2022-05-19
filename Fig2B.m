% Figure 2 B, recapituale experimental results from Davis et. Al 2015
% Run simulations for MBNL1 overexpression in basal conditions and MBNL1
% knockdown in a high TGFB/AngII condition
%digitized data from Davis et. al 2015 is loaded from .mat files, acquired
%using Matlab's Digitize2 function
%Last update: 3-7-2022 by AN
%% Load initial experimental data from Davis et. Al 2015
clc;
clear;
modelfname='fibroblastMBNL1'; %declare model file name
%load digitized data from Davis et. Al 2015
calcData=load('calineurinDigitizedDataDavis.mat') %calcineurin data in MBNL1 OE
srfData=load('srfDigitizedDataDavis.mat') %SRF data in MBNL1 KO
%SRF Data
srfp=srfData.df([1,2],2); %means
srferr=srfData.df([3,4],2)-srfData.df([1,2],2); %standard error, subtract mean from upper error bar value
srferr=srferr*sqrt(3); %convert to standard dev, experimental N=3
srferr=srferr/max(srfp); %min-max normalization
srfp=srfp/max(srfp); %min-max normalization
%Caclineurin Data
cnap=calcData.df([2,1],2); %means
cnaerr=calcData.df([4,3],2)-calcData.df([2,1],2); %standard error, subtract mean from upper error bar value
cnaerr=cnaerr*sqrt(3); %convert to standard dev, experimental N=3
cnaerr=cnaerr/max(cnap); %min-max normalization
cnap=cnap/max(cnap); %min-max normalization

%% Start esemble modeling section
%if exist, delete old ODE file before generating new file
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
util.textwrite('ODEfun.m',commandLine); %write the ODE file to working dir
%set random seed, using random seed 0 for all sims in this study
randn('seed',0);
%pre-acclocate array for ensemble results
ensemble=[];
COV= 0.0331; %Input COV model
%start ensemble for loop
for k = [1:150]
%declare y values for knockout and overexpression
ko=0;   
oe=1;
%random sampling for raised inputs
high_input=normrnd(0.6,COV*0.6)
while (high_input < 0) | (high_input > 1) %resample if values less than 0, greater than 1
    high_input=normrnd(0.6,COV*0.6);
end  
%keep track of simulation number
disp(['simulation ',num2str(k)])
%declare simulation conditions
alts={'w(1)=high_input; w(2)=high_input;','w(1)=high_input; w(2)=high_input; ymax(114)=ko;','','ymax(114)=oe; y0(114)=oe; tau(114)=100000;'}
%random sampling for basal inputs
rand_inputs=normrnd(0.1,COV*0.1,[1,9])
while any(rand_inputs < 0) %resample if values less than 0
    rand_inputs=normrnd(0.1,COV*0.1,[1,9])
end
%random sampling for mechanical input
 rand_mech=normrnd(0.725,COV*0.725);
    while (rand_mech < 0) | (rand_mech > 1)%resample if values less than 0, greater than 1
    rand_mech=normrnd(0.725,COV*0.725);
    end  
%pre-allocate array to capture values for each simualted experiment
yEnd=[];
for i = 1:length(alts) %iterate over each experimental simualtion condition
    alt=alts{i}; %obtain alternate parameters
% set up simulation params
tspan = [0 500]; % run out to ss
options = [];
[w,n,EC50,tau,ymax,y0] = paramList{:};
%get random inputs for other model inputs, AngII/TGFB will be overwirrten
%if alts is called
w([1,2,4:10])=rand_inputs;
w(3)=rand_mech;
%alter params to simulate experimental conditions
eval(alt);
rpar=[w;n;EC50];
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yEnd(:,i) = real(y(end,:)'); %get equilibrium values for all nodes in simulation
end
%store results from iteration in total ensemble array for outputs of
%interest, 116=SRFmRNA, 47=Calcineurin
ensemble(k,:)=[yEnd(116,1),yEnd(116,2),yEnd(47,3),yEnd(47,4)];
end
%% Save the simulation results
save('ensembleFig2_2_2022.mat','ensemble')
%% Load sim results, allocate data for plotting, and min-max normalize data
load('ensembleFig2_2_2022.mat')
%simulation values
srfm=[ensemble(:,1),ensemble(:,2)];%get SRF sim values
srfm=srfm/max(mean(srfm)); %min-max normalization
cnam=[ensemble(:,3),ensemble(:,4)];% get calcineurin sim values
cnam=cnam/max(mean(cnam)); %min-max normalization
srfSEM = std(srfm);%standard dev
cnaSEM = std(cnam);%standard dev
%Generate figure
fig=figure;
%SRF subplot
subplot(2,1,1)
b=bar([srfp',mean(srfm)])
b.FaceColor=[29,145,192]/256
hold on
e1=errorbar([1,2],srfp',srferr)
e1.LineStyle='none';
e1.Color='black';
hold on
e=errorbar([3,4],mean(srfm),std(srfm))
e.LineStyle='none';
e.Color='black';
xlabel('Treatment')
ylabel('Normalized Expression')
title('SRFmRNA')
ylim([0,1.25])
set(gca, 'XTick', [1:4])
set(gca, 'XTickLabel', {'WT_E_x_p','MBNL1-KO_E_x_p','WT_M_o_d_e_l','MBNL1-KO_M_o_d_e_l'})
%Calcineurin subplot
subplot(2,1,2)
b=bar([cnap',mean(cnam)])
b.FaceColor=[29,145,192]/256
hold on
e1=errorbar([1,2],cnap',cnaerr)
e1.LineStyle='none';
e1.Color='black';
hold on
e=errorbar([3,4],mean(cnam),std(cnam))
e.LineStyle='none';
e.Color='black';
xlabel('Treatment')
ylabel('Normalized Expression')
title('Calcineurin')
ylim([0,1.7])
set(gca, 'XTick', [1:4])
set(gca, 'XTickLabel', {'WT_E_x_p','MBNL1-OE_E_x_p','WT_M_o_d_e_l','MBNL1-OE_M_o_d_e_l'})
saveas(fig,['Figure2B.png'])
saveas(fig,['Figure2B.svg'])
