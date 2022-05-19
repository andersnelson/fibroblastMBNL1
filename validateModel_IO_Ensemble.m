function [percentMatch, resultChart, BMatch] = validateModel_IO_Ensemble(modelfname, validationfname, inputW,COV)
% Calculates the percent agreement and generates resultChart variable 

%   inputs: 
%   - modelfname = model file name (.xlsx)
%   - validationfname = validation file name (.xlsx)
%   - inputW = mean value for basal inputs
%   - COV = input COV for random sampling
%   
%   outputs:
%   - percentMatch = percent agreement
%   - resultChart = chart containing results of each individual validation
%     simulation. 
%   - BMatch = boolean vector containing the result of each validation (1 =
%     correct, 0 = incorrect)

% delete the previously formed ODE to make sure it's rewritten
pwd = cd;
if exist([pwd '\ODEfun.m'],'file') == 2
    delete('ODEfun.m');
end

%set random seed
randn('seed',0);
% parse out model name (xls2Netflux needs it as an arg)
namepos = findstr('.xls', modelfname);
namestr = modelfname(1:namepos-1);
namestr = cellstr(namestr);
% generate ODE from model spreadsheet
[specID,reactionIDs,~,paramList,ODElist,~, error] = util.xls2Netflux(namestr,modelfname);
commandLine = util.exportODE2(specID,paramList,ODElist);
util.textwrite('ODEfun.m',commandLine); %generate ODE file
% set up simulation params
tspan = [0 500]; % run out to ss
options = [];
[w,n,EC50,tau,ymax,y0] = paramList{:};
% read the validation sheet
[~, txt, raw] = xlsread(validationfname);
% remove rows without data
noData = cellfun(@(x)isequal(x,'No Data'), txt(1:end, 7));
txt(noData, :) = []; 
noData = cellfun(@isempty, txt(1:end, 7));
txt(noData, :) = [];
assignin('base', 'txt', txt);
%index to columns of the validation sheet
input1 = txt(2:end, 2);
inputCode = txt(2:end, 3); %input code to be evaluated
measurement = txt(2:end,6); %output measured
outputSpec = txt(2:end, 4);
control = txt(2:end, 5);
validationIDs = txt(2:end, 1);
validationTags = txt(2:end, 8);

% convert species and rxn names to integer values to map name and reaction
% ID's to numerical integers, allowing the input code to be evaluated
% directly
for k = 1:length(specID)
    eval(strcat(specID{k},' = ',num2str(k),';'));
end
for i = 1:length(reactionIDs)
    eval(strcat(reactionIDs{i},' = ',num2str(i),';'));
end
for i = 1:length(validationIDs)
    if ~isempty(validationIDs{i}) 
        eval(strcat(validationIDs{i}, ' = ', num2str(i), ';'));
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

randis={};%random basal input values
randms={};%random mechanical input values
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
randis{q}=rand_input;
randms{q}=rand_mech1;
end

% loop over all validation simulations read in from the excel sheet
for i = 1:length(inputCode)
    ensemble=[];
    ensemble_fold=[];
    for q=1:150
     %acquire random inputs
    rand_inputs=randis{q};
    rand_mech=randms{q};
    %Control simulation
    [w,n,EC50,tau,ymax,y0] = paramList{:}; % reset params
    %assign normal-randomly sampled values to model input variables
    w([1,2,4:10])=rand_inputs;
    w(3)=rand_mech;
    disp(['testing relationship ', num2str(i),' # of ',num2str(length(inputCode))])
    %repack parameters
    rpar = [w;n;EC50];
    params = {rpar,tau,ymax,specID};
    [~,y] = ode15s(@ODEfun, tspan, y0, options, params);
    yStart = y(end,:)'; % use the "no input" steady state as control 
    
    %Experimental simulation
    [w,n,EC50,tau,ymax,y0] = paramList{:}; % reset params
    %assign normal-randomly sampled values to model input variables
    w([1,2,4:10])=rand_inputs;
    w(3)=rand_mech;
    %evaulate new change for experimental sim
    eval(inputCode{i});
    %repack params
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


    %Ouput activity change
    ensemble(q) = real(yEndL{i}(outputSpeciesIndex(i)))-real(yStartL{i}(outputSpeciesIndex(i)));
    %output fold change
    ensemble_fold(q)=real(yEndL{i}(outputSpeciesIndex(i)))/real(yStartL{i}(outputSpeciesIndex(i)));

    end
    activityChange=mean(ensemble)
    fold_change{i}=mean(ensemble_fold)
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

%% Output the results chart
resultChart = {validationIDs, input1, outputSpec, measurement, prediction', predChange', match', validationTags,fold_change'}; %create a cell array showing input, output, and whether they matched validation
header = {'ID', 'input' , 'output', 'measurement', 'prediction', 'predicted change', 'match', 'tag','fold_change'};
resultChart = horzcat(resultChart{:});
resultChart = vertcat(header, resultChart);
percentMatch = numMatching/length(measurement)*100;
disp(['Percent match is ',num2str(percentMatch),' with ',num2str(sum(BMatch)),'/',num2str(length(BMatch)),' matching.'])
delete('ODEfun.m');