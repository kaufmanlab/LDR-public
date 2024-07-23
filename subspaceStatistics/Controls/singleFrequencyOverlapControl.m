function output = singleFrequencyOverlapControl
%% OVERVIEW

% This function controls for the fact that variation in rotational plane
% may be due to estimation noise of a single plane, rather than true
% variation. To control for this, single trials of each condition were
% partitioned into halves, the rotational planes for either half estimated,
% and then overlap estimated. This is then compared to the values found by
% the actual datasets. 

%% Control for varation.

load('singleFrequencyOverlap');
load('ShenoyMonkeyDataSingleTrial');
load('TBFResults');

% For each monkey perform analysis.
for monkey = 1:size(ShenoyMonkeyDataSingleTrial,2)
    output(monkey).M1 = analyzeData( ...
        ShenoyMonkeyDataSingleTrial(monkey).M1, ...
        TBFResults(monkey).M1,singleFrequencyOverlap(monkey).M1);
    output(monkey).PMd = analyzeData( ...
        ShenoyMonkeyDataSingleTrial(monkey).PMd, ...
        TBFResults(monkey).PMd,singleFrequencyOverlap(monkey).PMd);
end

end

%% SUBFUNCTION FOR ANALYZING DATA
function out = analyzeData(singleTrial,results,overlap)
% Prune trials.
singleTrials = preprocessTrials(pruneRepeatsTrial(singleTrial));
% Loop through each basis fxn.
% Get the pairs of basis functions (for planes). 
IDs = pairID(results.params.param(:,1));
uniqueIDs = unique(IDs);
% Assign the number of folds.
k = 200;
% Loop through pairs to assign array.
for ID = uniqueIDs
    out(ID).estimates = zeros(size(results.loadings,2),k);
    for cond = 1:size(results.loadings,2)
        out(ID).useTensor(cond).tensor = zeros(size(results.loadings(1).matrix,1), ...
            length(find(IDs == ID)),k*2);
    end
end
% For k number of folds, estimate the magnitude of estimation error for
% a single space.
for cv = 1:k
    % Partition trials.
    [part1,part2] = cvConds(singleTrials,0.5);
    % Fit loading matrices.
    for cond = 1:size(results.loadings,2)
        loadings1(cond).matrix = part1(cond).matrix*pinv(results.params.rawBasisFxns.');
        loadings2(cond).matrix = part2(cond).matrix*pinv(results.params.rawBasisFxns.');
    end
    % Loop through responses.
    for ID = uniqueIDs
        % Get the usable indices.
        inds = find(IDs == ID);
        % For each condition, find the error.
        for cond = 1:size(results.loadings,2)
            out(ID).estimates(cond,cv) = ...
                getAlignmentIndex( ...
                loadings1(cond).matrix(:,inds)*results.basisFxns(:,inds).', ...
                loadings2(cond).matrix(:,inds));
            % Assign the subspaces to the tensor for later SEAs.
            out(ID).useTensor(cond).tensor(:,:,1+(cv-1)*2) = ...
                loadings1(cond).matrix(:,inds);
            out(ID).useTensor(cond).tensor(:,:,2+(cv-1)*2) = ...
                loadings2(cond).matrix(:,inds);
        end
    end
end
% Loop through pairs to quantify statistics.
for ID = uniqueIDs
    out(ID).mean = mean(out(ID).estimates(:));
    out(ID).std = std(out(ID).estimates(:));
    out(ID).pVal = ranksum(out(ID).estimates(:), ...
        overlap(ID).alignmentIndexForStats);
    % Calculate ROC-AUC.
    [out(ID).XCoord,out(ID).YCoord,~,~] = ...
        perfcurve([out(ID).estimates(:)*0; overlap(ID).alignmentIndexForStats*0+1], ...
        [out(ID).estimates(:); overlap(ID).alignmentIndexForStats],1);
    for coord = 1:length(out(ID).YCoord)
        if out(ID).YCoord(coord) < out(ID).XCoord(coord)
            storeXCoord = out(ID).XCoord(coord);
            storeYCoord = out(ID).YCoord(coord);
            out(ID).XCoord(coord) = storeYCoord;
            out(ID).YCoord(coord) = storeXCoord;
        end
    end
    out(ID).AUC = sum((out(ID).XCoord(2:end) ...
        - out(ID).XCoord(1:end-1)).*out(ID).YCoord(2:end));    
    % Quantify SEAs.
    out(ID).SEA = zeros(25,size(results.loadings,2));
    for cond = 1:size(results.loadings,2)
        out(ID).SEA(:,cond) = ...
            subspaceExcursionAngles(out(ID).useTensor(cond).tensor,25)*(360/(2*pi));
    end
    % Get the pVals.
    out(ID).pVals = zeros(25,1);
    for pVal = 1:25
        out(ID).pVals(pVal) = ...
            sum(out(ID).SEA(pVal,:) > overlap(ID).SEA(pVal))/size(out(ID).SEA,2);
    end
end
end

