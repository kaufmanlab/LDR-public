function output = getCorrOverlapWithjPCA
%% OVERVIEW

% This function correlates the overlap between rotational activity and the
% jPCA plane with jPCA variance. 

%% Find correlations.

% Load the data.
load('TBFResults');
load('jPCAResults');

% Loop through and find the planes.
for monkey = 1:size(TBFResults,2)
    [output(monkey).M1.corrCoeffs,output(monkey).M1.pVals, ...
    output(monkey).M1.jPCAVar,output(monkey).M1.overlap] ...
        = corrOverlapWithjPCA(TBFResults(monkey).M1.loadings, ...
        TBFResults(monkey).M1.basisFxns,TBFResults(monkey).M1.params.param, ...
        jPCAResults(monkey).M1.rotations);
    [output(monkey).PMd.corrCoeffs,output(monkey).PMd.pVals, ...
    output(monkey).PMd.jPCAVar,output(monkey).PMd.overlap] ...
        = corrOverlapWithjPCA(TBFResults(monkey).PMd.loadings, ...
        TBFResults(monkey).PMd.basisFxns,TBFResults(monkey).PMd.params.param, ...
        jPCAResults(monkey).PMd.rotations);
end

end