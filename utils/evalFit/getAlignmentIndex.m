function alignmentIndex = getAlignmentIndex(data,subspace)
%% OVERVIEW

% This function takes data and a subspace and returns the subspace overlap
% of the data with the subspace, quantifying how much of the data lives in
% that subspace.

%% Quantify the variance explained.

% Remove parts of the data not approximated, assuming that approximation is
% uniform across neurons but not time.
inds = find(isnan(data(1,:)));
data(:,inds) = [];

% Quantify the alignment index.
alignmentIndex = 1 -sum(var(data - subspace*pinv(subspace)*data,0,2))/ ...
    sum(var(data,0,2));
alignmentIndex(alignmentIndex < 0) = 0;


end

