function output = TBFmouse
%% OVERVIEW

% This function performs Temporal Basis Factorization to identify temporal 
% basis functions that are conserved across conditions, identifying the 
% subspaces in which they exist.

%% Perform temporal basis factorization.

% Load the dataset(s). 
load('beginToEndTransport22');

% Assign the parameter of the basis function number (I'll figure out how to
% cross-validate this later).
basisFxnNum = 8;

% Perform temporal basis factorization.
[output.basisFxns,output.loadings] = getBasisFxns(data.M1,basisFxnNum);

% Add a field to quantify the variance explained.
output.varExplained = ...
    getVarExplained(data.M1,useFactorization(output.loadings,output.basisFxns),'ind');
output.recoverability = getRecoverability( ...
        data.M1,output.loadings,output.basisFxns);

end