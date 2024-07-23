function output = analyzeDimDiffs
%% OVERVIEW

% This function analyzes whether differences in identified dimensionalities
% have any connection to reach direction. This is done by, for each reach
% direction (left-right), calculating the likelihood of getting the other
% direction's dimensionalities from the distribution of the chosen side. 

%% Anaylze differences.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('dimCVResults');

% For each dataset, analyze.
for monkey = 1:size(ShenoyMonkeyData,2)
    angles = extractAngles(pruneRepeats(ShenoyMonkeyData(monkey).Kinematics));
    output(monkey).M1 = compareMultinomials(dimCVResults(monkey).M1.maxDimCond, ...
        angles(1,:));
    output(monkey).PMd = compareMultinomials(dimCVResults(monkey).PMd.maxDimCond, ...
        angles(1,:));
end

end

%% COMPARE MULTINORMIALS

function output = compareMultinomials(dimResults,angles)
% Segrate left from right.
left = cos((2*pi/360)*angles) < -1/sqrt(2);
right = cos((2*pi/360)*angles) > 1/sqrt(2);
% Assign p values.
output.p = ranksum(dimResults(left),dimResults(right));
output.leftMean = mean(dimResults(left));
output.rightMean = mean(dimResults(right));
end