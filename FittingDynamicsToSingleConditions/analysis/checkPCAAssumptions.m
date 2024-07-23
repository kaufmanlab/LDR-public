function output = checkPCAAssumptions
%% OVERVIEW

% This function takes in the distribution of subspaces found by
% condition-specific dynamics, and checks whether they are more more
% orthogonal than would be expected by random, low-D realizations of a
% subspace. That is, we are testing against the hypothesis that the reason
% we get dissimilar subspaces is that we are just limitted by a noise
% floor, and that each subspace is drawn from the eigenvectors of the same
% covariance matrix, but just that the space is high dimensional enough
% that there is dissimilarity. 

%% Check the assumptions of PCA.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);

% Load the results.
load('condSpecificDynamics');

% For each dataset check the assumption.
for monkey = 1:size(ShenoyMonkeyData,2)
    output(monkey).M1 = analysis(ShenoyMonkeyData(monkey).M1, ...
        condSpecificDynamics(monkey).M1.alignmentIndex, ...
        condSpecificDynamics(monkey).M1.subspaces);
    output(monkey).PMd = analysis(ShenoyMonkeyData(monkey).PMd, ...
        condSpecificDynamics(monkey).PMd.alignmentIndex, ...
        condSpecificDynamics(monkey).PMd.subspaces);
end

end

%% SUBFUNCTION FOR ANALYZING DATA.

function output = analysis(data,alignmentIndices,subspaces)
data = pruneRepeats(data);
% Assign the original alignment indicies.
output.alignmentIndices = alignmentIndices;
output.meanAlignmentIndex = mean(alignmentIndices(:));
output.stdAlignmentIndex = std(alignmentIndices(:));
% Generate a null distribution with the loading method.
nullIndices1 = PCANullOverlap(data,subspaces,'loading');
output.nullIndices1 = nullIndices1;
output.meanNullIndex1 = mean(nullIndices1(:));
output.stdNullIndex1 = std(nullIndices1(:));
% Check the pval.
output.pval1 = ranksum(nullIndices1(:),alignmentIndices(:));
% Generate a null distribution with the randDraw method.
nullIndices2 = PCANullOverlap(data,subspaces,'randDraw');
output.nullIndices2 = nullIndices2;
output.meanNullIndex2 = mean(nullIndices2(:));
output.stdNullIndex2 = std(nullIndices2(:));
% Check the pval.
output.pval2 = ranksum(nullIndices2(:),alignmentIndices(:));
end