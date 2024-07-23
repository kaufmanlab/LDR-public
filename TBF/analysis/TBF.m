function output = TBF
%% OVERVIEW

% This function performs Temporal Basis Factorization to identify the
% conserved eigenvalues of motor cortex dynamics during movement. This is
% done by identifying temporal basis functions that are conserved across
% conditions, identifying the subspaces in which they exist, and then
% fitting parametric functions that describe the temporal basis
% function/the dynamics they describe. 

%% Perform temporal basis factorization.

% Load the datasets. 
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);

% Load the number of basis functions.
load('basisFxnNum');

% Perform temporal basis factorization.
for monkey = 1:size(ShenoyMonkeyData,2)
    [output(monkey).M1.basisFxns,output(monkey).M1.loadings, ...
        output(monkey).M1.params] = eigTransform( ...
        pruneRepeats(ShenoyMonkeyData(monkey).M1), ...
        basisFxnNum(monkey).M1.maxDim);
    [output(monkey).PMd.basisFxns,output(monkey).PMd.loadings, ...
        output(monkey).PMd.params] = eigTransform( ...
        pruneRepeats(ShenoyMonkeyData(monkey).PMd), ...
        basisFxnNum(monkey).PMd.maxDim);
end

end