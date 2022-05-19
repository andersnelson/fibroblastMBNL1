% Figures 4 and 5, Simulate overexpressions of MBNL1 targets in a high
% MBNL1 signaling context and measure expressions for nodes of interest.
% Map cell signaling from MBNL1 target to aSMA
%Last update: 3-8-2022 by AN
clc;
clear;
%% model setup
% delete the previously formed ODE if it exists and generate new temp model 
% ODE file
if exist([pwd '\ODEfun.m'],'file') == 2
delete('ODEfun.m');
end
%declare model file
modelfname='fibroblastMBNL1.xlsx';
% parse out model name (xls2Netflux needs it as an arg)
namepos = findstr('.xls', modelfname);
namestr = modelfname(1:namepos-1);
namestr = cellstr(namestr);
% generate ODE from model spreadsheet
[specID,reactionIDs,~,paramList,ODElist,~, error] = util.xls2Netflux(namestr,modelfname);
commandLine = util.exportODE2(specID,paramList,ODElist);
util.textwrite('ODEfun.m',commandLine);

%% Run simulations for Figure 5
%nodes to plot
plottedNodes1=[115,85,48,87]; %nodes to plot in MK2 OE: SRF,NFAT,aSMA
plottedNodes2=[48,85,15,46,47,87]; %nodes to plot in NFAT OE: SRF,TRPC,Ca,calcineruin,aSMA
%run simulatons
tspan = [0 500]; % run out to steady-state
options = []; %options var for ODE solver
ensemble=[]; %allocate array for ensemble simulation results
COV=0.0331; %input COV
randn('seed',0); %set random seed
%% simualtion 1, MK2 OE
%initiate ensemble for loop
for h=1:150
%Control simulation
%random inputs
rand_inputs=normrnd(0.1,COV*0.1,[1,7]) 
while any(rand_inputs < 0 | 1 < rand_inputs) %resample if out of range
    rand_inputs=normrnd(0.1,COV*0.1,[1,7])
end    
rand_mech=normrnd(0.725,COV*0.725)
while rand_mech < 0 | 1 < rand_mech %resample if out of range
    rand_mech=normrnd(0.725,COV*0.725)
end  
rand_raised=normrnd(0.25,COV*0.25,[1,2])
while any(rand_raised < 0 | 1 < rand_raised) %resample if out of range
    rand_raised=normrnd(0.25,COV*0.25,[1,2])
end    
%unpack params
[w,n,EC50,tau,ymax,y0] = paramList{:};
%set randomly sampled input values
w([1,2])=rand_raised; %raised TGFB and AngII
w([4:10])=rand_inputs; %random basal input
w(3)=rand_mech;%random mech input
rpar=[w;n;EC50];
%high MBNL1
y0(114)=1;
tau(114)=100000;
%pack params
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yControl = real(y(end,:)'); 
%Experimental sim
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2])=rand_raised;%raised TGFB and AngII
w([4:10])=rand_inputs;
w(3)=rand_mech;
rpar=[w;n;EC50];
%MBNL1 and MK2 OE
y0(114)=1;
tau(114)=100000;
y0(115)=1;
tau(115)=100000;
%re-run exp. sim with new params
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yExp = real(y(end,:)'); 
delta=yExp(plottedNodes1)-yControl(plottedNodes1);
ensemble(h,:)=delta;
end
%% Generate Figure
values=mean(ensemble);
errors=(std(ensemble));
fig=figure;
hb=bar(values)
xticklabels({'MK2','SRF','NFAT','\alphaSMA'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14)
title('MK2 Overexpression','FontSize',18)
ylabel('Normalized Expression Increase','FontSize',16)
xtickangle(30)
hb.FaceColor = 'flat';
hb.CData(1,:) = [8,29,88]/256;
hb.CData(2,:) = [37,52,148]/256;
hb.CData(3,:) = [34,94,168]/256;
hb.CData(4,:) = [237,248,177]/256;
hold on
e1=errorbar(values,errors)
e1.LineStyle='none';
e1.Color='black'
saveas(fig,['Fig5_MK2.pdf'])

%% simualtion 2, NFAT OE
tspan = [0 500]; % run out to ss
options = []; %arg for ODE sovler
ensemble=[]; %allocate array for ensemble sim. results
COV=0.0331; %input COV
randn('seed',0); %set random seed
%initiate ensemble sim. for loop
for h=1:150
%control sim
rand_inputs=normrnd(0.1,COV*0.1,[1,7])
while any(rand_inputs < 0 | 1 < rand_inputs)
    rand_inputs=normrnd(0.1,COV*0.1,[1,7]) %resample if out of range
end    
rand_mech=normrnd(0.725,COV*0.725)
while rand_mech < 0 | 1 < rand_mech
    rand_mech=normrnd(0.725,COV*0.725) %resample if out of range
end  
rand_raised=normrnd(0.25,COV*0.25,[1,2])
while any(rand_raised < 0 | 1 < rand_raised)
    rand_raised=normrnd(0.25,COV*0.25,[1,2]) %resample if out of range
end    
%unpack params
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2])=rand_raised; %raised TGFB and AngII
w([4:10])=rand_inputs;
w(3)=rand_mech;
rpar=[w;n;EC50];
%high MBNL1
y0(114)=1;
tau(114)=100000;
%repack params
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yControl = real(y(end,:)'); 
%Experimental sim
%unpack params
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2])=rand_raised; %raised AngII and TGFB
w([4:10])=rand_inputs; %raised bsal inputs
w(3)=rand_mech; %raised mech input
rpar=[w;n;EC50];
%MBNL1 and NFAT OE
y0(114)=1;
tau(114)=100000;
y0(48)=1;
tau(48)=100000;
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yExp = real(y(end,:)'); 
delta=yExp(plottedNodes2)-yControl(plottedNodes2);
ensemble(h,:)=delta;
end
%% Generate Figure
values=mean(ensemble);
errors=(std(ensemble));
fig=figure;
hb=bar(values)
xticklabels({'NFAT','SRF','TRPC','Calcium','calcineruin','\alphaSMA'})
title('NFAT Overexpression','FontSize',18)
ylabel('Normalized Expression Increase','FontSize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14)
xtickangle(30)
hb.FaceColor = 'flat';
hb.CData(1,:) = [34,94,168]/256;
hb.CData(2,:) = [37,52,148]/256;
hb.CData(3,:) = [29,145,192]/256;
hb.CData(4,:) = [65,182,196]/256;
hb.CData(5,:) = [127,205,187]/256;
hb.CData(6,:) = [237,248,177]/256;
hold on
e1=errorbar(values,errors)
e1.LineStyle='none';
e1.Color='black'
saveas(fig,['Fig5_NFAT.pdf'])

%% Figure 4, Hippo, Sox9, Runx
% model setup, clear old temp. ODE file
if exist([pwd '\ODEfun.m'],'file') == 2
delete('ODEfun.m');
end
%declare model file
modelfname='fibroblastMBNL1.xlsx';
% parse out model name (xls2Netflux needs it as an arg)
namepos = findstr('.xls', modelfname);
namestr = modelfname(1:namepos-1);
namestr = cellstr(namestr);
% generate ODE from model spreadsheet
[specID,reactionIDs,~,paramList,ODElist,~, error] = util.xls2Netflux(namestr,modelfname);
commandLine = util.exportODE2(specID,paramList,ODElist);
util.textwrite('ODEfun.m',commandLine);

%% Run simulations for figure 4
%nodes to plot
plottedNodes1=[117,61,62,120,87]; %nodes to plot in Sox9OE, Sox9, PI3k, Akt, Foxo3,aSMA
plottedNodes2=[119,29,87]; %nodes to plot in Runx and Cbfb OE,Runx1, Cbfb, smad3, aSMA 
plottedNodes3=[113,29,30,87];%nyap1 OE, nYap, Smad3, smad7, aSMA
%simulatons
tspan = [0 500]; % run out to ss
options = []; %options var for ODE solver

%% simualtion 1, Sox9 OE
ensemble=[]; %preallocate array for ensemble sim. values
COV=0.0331; %input COV
randn('seed',0); %set random seed
%initiate ensemble sim for loop
for h=1:150
%Control sim
rand_inputs=normrnd(0.1,COV*0.1,[1,7])
while any(rand_inputs < 0 | 1 < rand_inputs) %resample if out of range
    rand_inputs=normrnd(0.1,COV*0.1,[1,7])
end    
rand_mech=normrnd(0.725,COV*0.725)
while rand_mech < 0 | 1 < rand_mech %resample if out of range
    rand_mech=normrnd(0.725,COV*0.725)
end  
rand_raised=normrnd(0.25,COV*0.25,[1,2])
while any(rand_raised < 0 | 1 < rand_raised) %resample if out of range
    rand_raised=normrnd(0.25,COV*0.25,[1,2])
end    
%Control sim
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2])=rand_raised; %raised AngII and TGFB
w([4:10])=rand_inputs;
w(3)=rand_mech;
rpar=[w;n;EC50];
%high MBNL1
y0(114)=1;
tau(114)=100000;
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yControl = real(y(end,:)'); 
%Experimental sim
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2])=rand_raised; %raised AngII and TGFB
w([4:10])=rand_inputs;
w(3)=rand_mech;
rpar=[w;n;EC50];
%High MBNL1 and Sox9 OE
y0(114)=1;
tau(114)=100000;
y0(117)=1;
tau(117)=100000;
%repack params
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yExp = real(y(end,:)'); 
delta=yExp(plottedNodes1)-yControl(plottedNodes1)
ensemble(h,:)=delta;
end
%% Generate Figure
values=mean(ensemble);
errors=(std(ensemble));
fig=figure;
hb=bar(values)
xticklabels({'Sox9','PI3K','Akt','Foxo3','\alphaSMA'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14)
title('Sox9 Overexpression','FontSize',18)
ylabel('Normalized Expression Increase','FontSize',16)
xtickangle(30)
hb.FaceColor = 'flat';
hb.CData(1,:) = [199,233,180]/256;
hb.CData(2,:) = [127,205,187]/256;
hb.CData(3,:) = [65,182,196]/256;
hb.CData(4,:) = [29,145,192]/256;
hb.CData(5,:) = [237,248,177]/256;
hold on
e1=errorbar(values,errors)
e1.LineStyle='none';
e1.Color='black'
saveas(fig,['Fig4_Sox9.pdf'])

%% simualtion 2, Runx1
ensemble=[]; %preallocate array for ensemble sim. values
COV=0.0331; %input COV
randn('seed',0); %set random seed
for h=1:150
%Control sim
rand_inputs=normrnd(0.1,COV*0.1,[1,7]) 
while any(rand_inputs < 0 | 1 < rand_inputs) %resample if out of range
    rand_inputs=normrnd(0.1,COV*0.1,[1,7])
end    
rand_mech=normrnd(0.725,COV*0.725)
while rand_mech < 0 | 1 < rand_mech %resample if out of range
    rand_mech=normrnd(0.725,COV*0.725)
end  
rand_raised=normrnd(0.25,COV*0.25,[1,2])
while any(rand_raised < 0 | 1 < rand_raised) %resample if out of range
    rand_raised=normrnd(0.25,COV*0.25,[1,2])
end    
%Control sim
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2])=rand_raised;%raised AngII and TGFB
w([4:10])=rand_inputs;
w(3)=rand_mech;
rpar=[w;n;EC50];
%high MBNL1
y0(114)=1;
tau(114)=100000;
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yControl = real(y(end,:)'); 
%Experimental sim
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2])=rand_raised; %raised AngII and TGFB
w([4:10])=rand_inputs;
w(3)=rand_mech;
rpar=[w;n;EC50];
%High MBNL1 and Runx1 OE
y0(114)=1;
tau(114)=100000;
y0([119])=1;
tau([119])=100000;
%repack params
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yExp = real(y(end,:)'); 
delta=yExp(plottedNodes2)-yControl(plottedNodes2)
ensemble(h,:)=delta;
end

%% Generate Figure
values=mean(ensemble);
errors=(std(ensemble));
fig=figure;
hb=bar(values)
xticklabels({'Runx1','Smad3','\alphaSMA'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14)
title('Runx1 Overexpression','FontSize',18)
ylabel('Normalized Expression Increase','FontSize',16)
xtickangle(30)
hb.FaceColor = 'flat';
hb.CData(1,:) = [127,205,187]/256;
hb.CData(2,:) = [65,182,196]/256;
hb.CData(3,:) = [237,248,177]/256;
hold on
e1=errorbar(values,errors)
e1.LineStyle='none';
e1.Color='black'
saveas(fig,['Fig4_Runx.pdf'])



%% simualtion 3, Yap1
ensemble=[]; %allocate array for ensemble sim values
COV=0.0331; %input COV
randn('seed',0); %set random seed
%initiate ensemble sim. for loop
for h=1:150
%Control sim
rand_inputs=normrnd(0.1,COV*0.1,[1,7])
while any(rand_inputs < 0 | 1 < rand_inputs) %resample if out of range
    rand_inputs=normrnd(0.1,COV*0.1,[1,7])
end    
rand_mech=normrnd(0.725,COV*0.725)
while rand_mech < 0 | 1 < rand_mech %resample if out of range
    rand_mech=normrnd(0.725,COV*0.725)
end  
rand_raised=normrnd(0.25,COV*0.25,[1,2])
while any(rand_raised < 0 | 1 < rand_raised) %resample if out of range
    rand_raised=normrnd(0.25,COV*0.25,[1,2])
end    
%Control sim
%upack params
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2])=rand_raised; %raised AngII and TGFB
w([4:10])=rand_inputs;
w(3)=rand_mech;
rpar=[w;n;EC50];
%high MBNL1
y0(114)=1;
tau(114)=100000;
%repack params
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yControl = real(y(end,:)'); 
%Experimental sim
%upack params
[w,n,EC50,tau,ymax,y0] = paramList{:};
w([1,2])=rand_raised; %raised AngII and TGFB
w([4:10])=rand_inputs;
w(3)=rand_mech;
rpar=[w;n;EC50];
%High MBNL1 and Yap1 OE
y0(114)=1;
tau(114)=100000;
y0(113)=1;
tau(113)=100000;
%repack params
params = {rpar,tau,ymax,specID};
%run simualtion
[~,y] = ode15s(@ODEfun, tspan, y0, options, params);
yExp = real(y(end,:)'); 
delta=yExp(plottedNodes3)-yControl(plottedNodes3)
ensemble(h,:)=delta;
end
%% Generate figure
values=mean(ensemble);
errors=(std(ensemble));
fig=figure;
hb=bar(values)
xticklabels({'nYap1', 'Smad3', 'smad7', 'aSMA'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14)
title('Nuclear Yap1 Overexpression','FontSize',18)
ylabel('Normalized Expression Increase','FontSize',16)
xtickangle(30)
hb.FaceColor = 'flat';
hb.CData(1,:) = [127,205,187]/256;
hb.CData(2,:) = [29,145,192]/256;
hb.CData(3,:) = [34,94,168]/256;
hb.CData(4,:) = [237,248,177]/256;
hold on
e1=errorbar(values,errors);
e1.LineStyle='none';
e1.Color='black'
saveas(fig,['Fig4_Yap.pdf'])















