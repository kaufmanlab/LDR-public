function output = profileDynamics
%% OVERVIEW

% This function explores the dynamics on each condition of a dataset by
% first fitting dynamics to each condition independently, and quantifying
% the distribution of variance explained in each condition. Afterwards,
% three metrics are quantified. First, the variance in eigenvalues is
% quantified. Second, the overlap between subspaces is quantified. Third,
% the variance explained by substituting condition-specific eigenvalues
% with the average eigenvalues is quantified. 

%% Profile dynamics.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);

% Load the cross-validation statistics.
load('dimCVResults');

% For each dataset fit dynamics.
for monkey = 1:size(ShenoyMonkeyData,2)
    output(monkey).M1 = analysis(ShenoyMonkeyData(monkey).M1, ...
        dimCVResults(monkey).M1.maxDim);
    output(monkey).PMd = analysis(ShenoyMonkeyData(monkey).PMd, ...
        dimCVResults(monkey).PMd.maxDim);
end

end

%% FUNCTION FOR COLLECTIONG STATS.

function output = analysis(data,dimUse)

% Prune repeats.
data = pruneRepeats(data);

% Fit condition-specific dynamics.
[subspaces,eigVals,params,~,sim,basisFxns] = idDynamics(data,dimUse);
output.eigVals = eigVals;
output.subspaces = subspaces;
output.averageEigVals = averageEigVals(eigVals,'mean');
output.sim = sim;
output.basisFxns = basisFxns;
output.params = params;

% Quantify the variance explained.
output.VE = getVarExplained(data,sim,'ind');

% Quantify the variance in the eigenvalues.
output.eigValVar = var(eigVals,0,2);
output.eigValDist = eigValDistance(eigVals,0);
output.paramVar = var(params.param,0,3);

% Quantify the subspace overlap.
output.alignmentIndex = zeros(size(data,2));
for cond1 = 1:size(data,2)
    for cond2 = 1:size(data,2)
        output.alignmentIndex(cond1,cond2) ...
            = getAlignmentIndex(data(cond1).matrix,subspaces(cond2).matrix);
    end
end

% Quantify the temporal overlap.
output.temporalSim = zeros(size(data,2));
basisFxns = basisFxns(20:end,:,:);
for cond1 = 1:size(data,2)
    for cond2 = 1:size(data,2)
        output.temporalSim(cond1,cond2) ...
            = 1 - sum(var(squeeze(basisFxns(:,:,cond1)) ...
            *pinv(squeeze(basisFxns(:,:,cond1))) ...
            *squeeze(basisFxns(:,:,cond2)) - squeeze(basisFxns(:,:,cond2)),0,2))/ ...
            sum(var(squeeze(basisFxns(:,:,cond2)),0,2));
    end
end

end