function tangling = findTanglingAngle(input,others)
%% OVERVIEW

% This function, for each input point, finds the closet point to that point
% then measures the cosine angle between the derivatives at those points.

%% Find tangling.

% Loop over time points.
tangling = zeros(size(input,1)-1,1);
for timepoint = 1:size(input,1)-1
    % Identify the current point.
    current = input(timepoint,:);
    [minD,closest] = min(sum((current-others(1:end-1,:)).^2,2).^0.5);
    % Get the derivatives.
    currentDiff = input(timepoint+1,:) - input(timepoint,:);
    others = diff(others);
    others(timepoint,:) = [];
    tangling(timepoint) = currentDiff*others(closest,:).' ...
        /((norm(currentDiff,2)*norm(others(closest,:),2)));
end
tangling(tangling < 0) = 0;

end