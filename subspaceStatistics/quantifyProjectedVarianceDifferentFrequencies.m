function output = quantifyProjectedVarianceDifferentFrequencies(loadings,basisFxns,params)
%% OVERVIEW

% This function compares the subspaces occupied by different frequency
% responses on pairs of conditions by quantifying their alignment index,
% along with the smallest angle between the two subspaces. 

%% Compare subspaces.

% Get the pairs of basis functions (for planes). 
IDs = pairID(params(:,1));
uniqueIDs = unique(IDs);

% Loop through pairs.
for ID = uniqueIDs
    % Get the usable indices.
    inds = find(IDs == ID);
    % Loop through the IDs again.
    for compID = uniqueIDs
        compInds = find(IDs == compID);
        %% Alignment index.
        % Assign an array for the metric.
        output(ID).compare(compID).alignmentIndex = zeros(size(loadings,2),size(loadings,2));
        % Loop over all pairs of conditions.
        for cond1 = 1:size(loadings,2)
            % You have to loop naively as alignment index is not symmetric. 
            for cond2 = 1:size(loadings,2)
                output(ID).compare(compID).alignmentIndex(cond1,cond2) ...
                    = getAlignmentIndex( ...
                    loadings(cond1).matrix(:,inds)*basisFxns(:,inds).', ...
                    loadings(cond2).matrix(:,compInds));
            end
        end
        output(ID).compare(compID).meanAlignmentIndex = mean(output(ID).compare(compID).alignmentIndex(:));
        output(ID).compare(compID).stdAlignmentIndex = std(output(ID).compare(compID).alignmentIndex(:));
        output(ID).compare(compID).useDist = output(ID).compare(compID).alignmentIndex;
        if ID == compID
            for cond = 1:size(loadings,2)
                output(ID).compare(compID).useDist(cond,cond) = nan;
            end
        end
        output(ID).compare(compID).useDist = output(ID).compare(compID).useDist(:);
        output(ID).compare(compID).useDist(isnan(output(ID).compare(compID).useDist)) = [];
        %% Smallest angle.
        % Assign an array for the metric.
        output(ID).compare(compID).angle = zeros(size(loadings,2),size(loadings,2));
        % Loop over all pairs of conditions.
        for cond1 = 1:size(loadings,2)
            % You have to loop naively as alignment index is not symmetric. 
            for cond2 = 1:size(loadings,2)
                angles = subspacea( ...
                    loadings(cond1).matrix(:,inds), ...
                    loadings(cond2).matrix(:,compInds));
                output(ID).compare(compID).angle(cond1,cond2) = angles(1)*360/(2*pi);
            end
        end
        output(ID).compare(compID).meanAngle = mean(output(ID).compare(compID).angle(:));
        output(ID).compare(compID).stdAngle = std(output(ID).compare(compID).angle(:));
    end
end

% Loop through pairs.
for ID = uniqueIDs
    % Loop through the IDs again.
    for compID = uniqueIDs
        %% Alignment index.
        % Assign an array for the metric.
        output(ID).compare(compID).alignmentIndexPVal = ...
            ranksum(output(ID).compare(ID).useDist(:), ...
            output(ID).compare(compID).useDist(:));
        % Calculate ROC-AUC.
        [output(ID).compare(compID).XCoord,output(ID).compare(compID).YCoord,~,~] = ...
            perfcurve([output(ID).compare(ID).useDist(:)*0; output(ID).compare(compID).useDist(:)*0+1], ...
            [output(ID).compare(ID).useDist(:); output(ID).compare(compID).useDist(:)],1);
        for coord = 1:length(output(ID).compare(compID).YCoord)
            if output(ID).compare(compID).YCoord(coord) < output(ID).compare(compID).XCoord(coord)
                storeXCoord = output(ID).compare(compID).XCoord(coord);
                storeYCoord = output(ID).compare(compID).YCoord(coord);
                output(ID).compare(compID).XCoord(coord) = storeYCoord;
                output(ID).compare(compID).YCoord(coord) = storeXCoord;
            end
        end
        output(ID).compare(compID).AUC = sum((output(ID).compare(compID).XCoord(2:end) ...
            - output(ID).compare(compID).XCoord(1:end-1)).*output(ID).compare(compID).YCoord(2:end));  
        %% Smallest angle.
        % Assign an array for the metric.
        output(ID).compare(compID).anglePVal = signrank( ...
            output(ID).compare(ID).angle(:),output(ID).compare(compID).angle(:));
    end
end

end