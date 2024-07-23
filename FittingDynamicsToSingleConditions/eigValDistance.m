function distance = eigValDistance(eigVals,compVals)
%% OVERVIEW

% This is a function defined to return the distance between two complex
% vectors in order to compare eigenvalues.

%% Compute distance.

if numel(compVals) == 1
    distance = zeros(size(eigVals,2));
    for cond1 = 1:size(eigVals,2)
        for cond2 = 1:size(eigVals,2)
            distance(cond1,cond2) = ((eigVals(:,cond1)-eigVals(:,cond2))' ...
                *(eigVals(:,cond1)-eigVals(:,cond2)))^0.5;
        end
    end
else
    distance = zeros(size(eigVals,2),1);
    for cond = 1:size(eigVals,2)
        distance(cond) = ((eigVals(:,cond)-compVals(:,cond))' ...
            *(eigVals(:,cond)-compVals(:,cond)))^0.5;
    end
end

