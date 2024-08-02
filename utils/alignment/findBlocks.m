function blocks = findBlocks(in)
%% OVERVIEW

% This function finds temporally-extended blocks in binary time series.

%% Find blocks.

% Assign the initial states.
cluster = 0;
inCluster = false;

% Loop over time periods.
for time = 1:length(in)
    if in(time) && ~inCluster
        cluster = cluster+1;
        inCluster = true;
        blocks(cluster).inds = time;
    elseif in(time) && inCluster
        blocks(cluster).inds = [blocks(cluster).inds time];
    elseif ~in(time) && inCluster
        inCluster = false;
    end
end

end

