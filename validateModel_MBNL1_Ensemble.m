function [percentMatch, resultChart1, inputArray] = validateModel_Int_Ensemble(modelfname,validationfname,COV)
% Calculate intermediate-outpu validations 
% First version edited 1-14-2020 by AN
% Last updated by AN 03-10-2022
%   inputs: 
%   - modelfname = model file name (.xlsx)
%   - validationfname = validation file name (.xlsx)
%   - alts = alterations in one or series of reactions, used in
%   iterative changes in the model 
%   
%   outputs:
%   - percentMatch = percent agreement
%   - resultChart = chart containing results of each individual validation
%     simulation. 
%   - BMatch = boolean vector containing the result of each validation (1 =
%     correct, 0 = incorrect)
%   - byClass = percent match by validation class

% delete the previously formed ODE to make sure it's rewritten
pwd = cd;
if exist([pwd '\ODEfun.m'],'file') == 2
    delete('ODEfun.m');
end
%set random seed
randn('seed',0);
% parse out model name (xls2Netflux needs it as an arg)
namepos = findstr('.xlsx', modelfname);
namestr = modelfname(1:namepos-1);
namestr = cellstr(namestr);

% generate ODE from model spreadsheet
[specID,reactionIDs,~,paramList,ODElist,~, error] = util.xls2Netflux(namestr,modelfname);
commandLine = util.exportODE2(specID,paramList,ODElist);
util.textwrite('ODEfun.m',commandLine);
% set up simulation params
tspan = [0 500]; % run out to steady-state
options = []; %options for ODE solver
[w,n,EC50,tau,ymax,y0] = paramList{:}; %get params
% read the validation sheet
[~, txt, raw] = xlsread(validationfname);
% remove rows without data
noData = cellfun(@(x)isequal(x,'No Data'), txt(1:end, 7));
txt(noData, :) = []; 
noData = cellfun(@isempty, txt(1:end, 7));
txt(noData, :) = [];
assignin('base', 'txt', txt);
subNetwork = txt(2:end, 1);
validationIDs = txt(2:end, 2);
input1 = txt(2:end, 3);
baseline = txt(2:end, 4); %these are baseline values, i.e. testing a KO in high TGFB context, the baseline is high tgfb
inputCode = txt(2:end, 5); % input codes to be evaluated for exp. sims
outputSpec=txt(2:end, 6);
measurement = txt(2:end,7); % measured output
description=txt(2:end,10);

%preallocate array for ensemble simulation results
ensembleSimResults=[];
% convert species and rxn names to integer values to map name and reaction
% ID's to numerical integers, allowing the input code to be evaluated
% directly
for k = 1:length(specID)
    eval([specID{k},' = ',num2str(k),';']);
end
for k = 1:length(reactionIDs)
    disp(k)
    eval([reactionIDs{k},' = ',num2str(k),';']);
    
end
for k = 1:length(validationIDs)
    if ~isempty(validationIDs{k}) 
        eval([validationIDs{k}, ' = ', num2str(k), ';']);
    end
end

% set validation threshold change
thresh = .001; % threshold, Ryall et al., 2012 set to 0.001 for sensitivity analysis
inc = {'Increase'};
dec = {'Decrease'};
noc = {'No Change'};

% determination of number of predictions matching references
numMatching = 0; % number of predictions consistent with the qualitative literature species behavior

% find indices of output species
outputSpeciesIndex = zeros(1, length(measurement));
for k = 1:length(outputSpec)
    [~,outputSpeciesIndex(k)] = ismember(outputSpec{k},specID);
end
inputArray=[];
randis={};
randms={};
%assign 150 random sets of simulation inputs
for q=1:150
    rand_input=normrnd(0.1,COV*0.1,[1,9]);
    while any(rand_input < 0 | 1 < rand_input)
    rand_input=normrnd(0.1,COV*0.1,[1,9]);
    end  
    rand_mech1=normrnd(0.725,COV*0.725);
    while (rand_mech1 < 0) | (rand_mech1 > 1)
    rand_mech1=normrnd(0.725,COV*0.725);
    end 
randis{q}=rand_input; %random basal inp
randms{q}=rand_mech1;
end

% loop over all validation simulations read in from the excel sheet
for i = 1:length(inputCode)
    ensemble=[];
    ensemble_fold=[];
    for q=1:150
    rand_inputs=randis{q};
    rand_mech=randms{q};
    
    %disp(['Simulation # ', num2str(i), ' of ',num2str(length(inputCode))]) % write the simulation number to the command line to track loop progress
    [w,n,EC50,tau,ymax,y0] = paramList{:}; % reset params
    w([1,2,4:10])=rand_inputs;
    w(3)=rand_mech;
    disp(['testing relationship ', num2str(i),' # of ',num2str(length(inputCode))])
    
%     % evaluate alternate parameters
%     eval(alts);
    if ~strcmp(baseline{i}, 'normal')
    eval(baseline{i});
    end
    rpar = [w;n;EC50];
    params = {rpar,tau,ymax,specID};
    [~,y] = ode15s(@ODEfun, tspan, y0, options, params);
    yStart = y(end,:)'; % use the "no input" steady state as control 
    
    
    % evaluate validation conditions from excel sheet
    [w,n,EC50,tau,ymax,y0] = paramList{:}; % reset params
    w([1,2,4:10])=rand_inputs;
    w(3)=rand_mech;
   
    eval(inputCode{i});
    if ~strcmp(baseline{i}, 'normal')
    eval(baseline{i});
    end
    rpar = [w;n;EC50];
    params = {rpar,tau,ymax,specID};
    [~,y] = ode15s(@ODEfun, tspan, y0, options, params);
    yEnd = y(end,:)'; % final species activation state
    % subtract final species activation from initial to determine the effect
    % the validation simulation has had on the species' activation with
    % respect to the defined threshold, then write out the qualitative
    % change of the species' activation state to the excel file
    yStartL{i} = yStart;
    yEndL{i} = yEnd;

    ensembleSimResults(q,i)=real(yEndL{i}(outputSpeciesIndex(i)))-real(yStartL{i}(outputSpeciesIndex(i)));
%     if isempty(control{i}) % if control validation defined
        ensemble(q) = real(yEndL{i}(outputSpeciesIndex(i)))-real(yStartL{i}(outputSpeciesIndex(i)));
        
        %output fold change as well as activity change
        fold_change{i}=real(yEndL{i}(outputSpeciesIndex(i)))/real(yStartL{i}(outputSpeciesIndex(i)));

inputArray(q,:)=rand_inputs;
    end
    activityChange=mean(ensemble);
    if activityChange > thresh % increase
        prediction{i} = 'Increase';
        predChange{i} = num2str(activityChange);
        fold_change{i}=num2str(fold_change{i});
        if isequal(inc,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1; %if the simulation matches the experimental validation put a 1 in the vector
        else
            match(i) = 0; %if the simulation does not match put a 0 in the matrix
        end
    elseif activityChange < -thresh % decrease
        prediction{i} = 'Decrease';
        predChange{i} = num2str(activityChange);
        if isequal(dec,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1;
        else
            match(i) = 0;
        end
    else % no change
        prediction{i} = 'No Change';
        predChange{i} = num2str(activityChange);
        if isequal(noc,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1;
        else
            match(i) = 0;
        end
    end
end

for j = 1:length(match)
    if match(j) == 1;
        match2{j} = 'yes';
        match3(j) = 1;
    else
        match2{j} = 'no';
        match3(j) = 0;
    end
end
match = match2;

BMatch = match3';

%% get total and percent validation for individual signaling subNetworks
%plot figure 2C
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
%generate figure 2C        
figure;
x=([1:length(usubNetwork)]);
values=[matched'./totals']*100;
b=bar(values)
b.FaceColor=[29,145,192]/256
set(gca,'XTickLabel',usubNetwork);
title('MBNL1 Model Experimental Valdiation')

%% Output ensemble simualation results
resultChart1 = {subNetwork,validationIDs, input1, baseline, outputSpec, measurement, prediction', predChange', match', fold_change',description}; %create a cell array showing input, output, and whether they matched validation
header = {'subNetwork','ID', 'input' ,'baseline' ,'output', 'measurement', 'prediction', 'predicted change', 'match', 'fold_change','description'};
resultChart1 = horzcat(resultChart1{:});
resultChart1 = vertcat(header, resultChart1);
percentMatch = numMatching/length(measurement)*100;
disp(['Percent match is ',num2str(percentMatch),' with ',num2str(sum(BMatch)),'/',num2str(length(BMatch)),' matching.'])
delete('ODEfun.m');

