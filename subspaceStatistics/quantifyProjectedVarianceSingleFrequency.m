function output = quantifyProjectedVarianceSingleFrequency(loadings,basisFxns,params)
%% OVERVIEW

% This functions compares the subspaces in which same frequency rotations
% live on different conditions. This is done by quantifying the projected
% variance overlap between the subspaces occupied on different conditions.

%% Quantify overlap.

% Get the pairs of basis functions (for planes). 
IDs = pairID(params(:,1));
uniqueIDs = unique(IDs);

% Loop through pairs.
for ID = uniqueIDs
    %% Alignment index.
    % Get the usable indices.
    inds = find(IDs == ID);
    % Assign an array for the metric.
    output(ID).alignmentIndex = zeros(size(loadings,2),size(loadings,2));
    % Loop over all pairs of conditions.
    for cond1 = 1:size(loadings,2)
        % You have to loop naively as alignment index is not symmetric. 
        for cond2 = 1:size(loadings,2)
            output(ID).alignmentIndex(cond1,cond2) = getAlignmentIndex( ...
                loadings(cond1).matrix(:,inds)*basisFxns(:,inds).', ...
                loadings(cond2).matrix(:,inds));
        end
    end
    % Create another array without the diagonals for statistical tests (the
    % first array can be regression purposes). 
    copyArray = output(ID).alignmentIndex;
    for cond = 1:size(loadings,2)
        copyArray(cond,cond) = nan;
    end
    copyArray = copyArray(:);
    copyArray(isnan(copyArray)) = [];
    output(ID).alignmentIndexForStats = copyArray;
    output(ID).meanAlignmentIndex = mean(copyArray);
    output(ID).stdAlignmentIndex = std(copyArray);
    % Get the excursion angles.
    useTensor = zeros(size(loadings(1).matrix,1), ...
        length(inds),size(loadings,2));
    for cond = 1:size(loadings,2)
        useTensor(:,:,cond) = loadings(cond).matrix(:,inds);
    end
    output(ID).SEA = subspaceExcursionAngles(useTensor,25)*(360/(2*pi));
    %% Smallest angle.
    % Get the usable indices.
    inds = find(IDs == ID);
    % Assign an array for the metric.
    output(ID).angle = zeros(size(loadings,2),size(loadings,2));
    % Loop over all pairs of conditions.
    for cond1 = 1:size(loadings,2)
        % You have to loop naively as alignment index is not symmetric. 
        for cond2 = 1:size(loadings,2)
            angles = subspacea( ...
                loadings(cond1).matrix(:,inds), ...
                loadings(cond2).matrix(:,inds));
            output(ID).angle(cond1,cond2) = angles(1)*360/(2*pi);
        end
    end
    % Create another array without the diagonals for statistical tests (the
    % first array can be regression purposes). 
    copyArray = output(ID).angle;
    for cond = 1:size(loadings,2)
        copyArray(cond,cond) = nan;
    end
    copyArray = copyArray(:);
    copyArray(isnan(copyArray)) = [];
    output(ID).angleForStats = copyArray;
    output(ID).meanAngle = mean(copyArray);
    output(ID).stdAngle = std(copyArray);
end

end