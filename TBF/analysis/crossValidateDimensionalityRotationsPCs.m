function output = crossValidateDimensionalityRotationsPCs
%% OVERVIEW

% This function cross-validates the dimensionality of each rotation
% independently. This is done by getting seperate estimates of each
% rotation from independent splits of the data, then the dimensionality of
% the subspace of each rotation is swept to maximize the variance explained
% in the alternate estimate of the rotation. As oppposed to the other
% similarly-named function, this function uses descriptive statistical
% versions of dimensionality, such as number of dimensions needed to
% capture 80% of the variance in the rotation, and the participation ratio.

%% Cross-validate.

% Load the dataset.
load('TBFResults');

% For each monkey, cross validate each region.
for monkey = 1:size(TBFResults,2)
    output(monkey).M1 = method(TBFResults(monkey).M1);
    output(monkey).PMd = method(TBFResults(monkey).PMd);
end

end

%% METHOD OF CROSSVALIDATING

function output = method(data)
% Get the frequencies present in the fxns.
uniqueFreqs = unique(data.params.param(:,1));
% Loop over rotations
for rot = 1:length(uniqueFreqs)
    % Get the indices of the rotation.
    findInds = find(data.params.param(:,1) == uniqueFreqs(rot));
    % Create an array for the rotations.
    rotationArray = zeros(size(data.loadings(1).matrix,1), ...
        length(find(~isnan(data.basisFxns(:,1))))*size(data.loadings,2));
    % Loop over conditions, populate array.
    for cond = 1:size(data.loadings,2)
        rotationArray(:,(1:length(find(~isnan(data.basisFxns(:,1))))) ...
            + (cond-1)*length(find(~isnan(data.basisFxns(:,1))))) ...
            = data.loadings(cond).matrix(:,findInds) ...
            * data.basisFxns(find(~isnan(data.basisFxns(:,1))),findInds).';
    end
    % Get the PCs of the data.
    [~,vars,~] = svd(rotationArray - mean(rotationArray,2),'econ');
    output(rot).vars = diag(vars)/sum(diag(vars));
    % Quantify the number of dimensions needed to explain 80% of the data.
    output(rot).dim80Percent = find(cumsum(diag(vars))/sum(diag(vars))>0.8);
    output(rot).dim80Percent = output(rot).dim80Percent(1);
    % Quantify the participation ratio.
    output(rot).PR = (sum(diag(vars)).^2)/sum(diag(vars).^2);
end
end

