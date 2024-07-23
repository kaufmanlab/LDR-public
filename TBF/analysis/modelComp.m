function output = modelComp
%% OVERVIEW

% This function just simply assembles and quantifies different model
% comparisons into one location, jPCA (same planes and frequencies), TBF
% (same frequencies, diff planes), single condition LDS (diff planes and
% frequencies, and compares a couple things about them. 

%% Perform temporal basis factorization.

% Load the datasets. 
load('TBFControls');
load('condSpecificDynamics');
load('jPCAResults');

% Perform temporal basis factorization.
for monkey = 1:size(TBFControls,2)
    output(monkey).M1 = analyzeStats(condSpecificDynamics(monkey).M1.VE.array, ...
        TBFControls(monkey).M1.data.varExplained, ...
        jPCAResults(monkey).M1.varExplained.array);
    output(monkey).PMd = analyzeStats(condSpecificDynamics(monkey).PMd.VE.array, ...
        TBFControls(monkey).PMd.data.varExplained, ...
        jPCAResults(monkey).PMd.varExplained.array);
end

end

%% Subfunction

function output = analyzeStats(condSpecific,tbf,jPCA)
output.condSpec = condSpecific;
output.tbf = tbf;
output.jPCA = jPCA;
end