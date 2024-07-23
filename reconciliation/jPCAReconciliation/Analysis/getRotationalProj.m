function output = getRotationalProj
%% OVERVIEW

% This function finds the projection of the data into the "jPCA rotational
% plane" using jPCA, along with the reconstruction due to rotations. The
% function then finds the variance in each rotational plane.

%% Find the rotational projection.

% Load the data.
load('TBFResults');

% Loop through and find the planes.
for monkey = 1:size(TBFResults,2)
    [output(monkey).M1.rotations, ...
        output(monkey).M1.reconstruction, ...
        output(monkey).M1.varExplained,output(monkey).M1.R2s] ...
        = getRotations(TBFResults(monkey).M1.loadings, ...
        TBFResults(monkey).M1.basisFxns,2);
    [output(monkey).PMd.rotations, ...
        output(monkey).PMd.reconstruction, ...
        output(monkey).PMd.varExplained,output(monkey).PMd.R2s] ...
        = getRotations(TBFResults(monkey).PMd.loadings, ...
        TBFResults(monkey).PMd.basisFxns,1);
end

end