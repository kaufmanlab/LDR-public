function output = fitSingleLDS
%% OVERVIEW

% This function finally puts reviewer 2 to rest.

%% Profile dynamics.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);

% Load the cross-validation statistics.
load('dimCVResults');

% For each dataset fit dynamics.
for monkey = 1:size(ShenoyMonkeyData,2)
    output(counter).M1 = analysis(ShenoyMonkeyData(monkey).M1);
    output(counter).PMd = analysis(ShenoyMonkeyData(monkey).PMd);
end

end

% This function fits a single LDS to neural activity. 

function out = analysis(neural)
% Copy neural data, truncate to peri-movement.
copyNeural1 = neural;
copyNeural2 = neural;
for cond = 1:size(neural,2)
    copyNeural1(cond).matrix = copyNeural1(cond).matrix(:,20:end-1);
    copyNeural2(cond).matrix = copyNeural2(cond).matrix(:,21:end);
end
% Fit an LDS. 
M = [copyNeural2().matrix]*pinv([copyNeural1().matrix]);
% Simulate each condition, calculate variance explained.
copyNeural = neural;
for cond = 1:size(neural,2)
    for t = 21:121
        copyNeural(cond).matrix(:,t) = M*copyNeural(cond).matrix(:,t-1);
    end
    out.varExplained(cond) = 1-sum(var(copyNeural(cond).matrix ...
        -neural(cond).matrix,0,2))/sum(var(neural(cond).matrix,0,2));
end
end