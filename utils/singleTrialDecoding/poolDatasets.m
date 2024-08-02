function pooledData = poolDatasets(d1,d2)
%% OVERVIEW

% This function is very basic, it just takes in two datasets and pools them
% to form a joint representation.

%% Pool

pooledData = d1;
for cond = 1:size(d1,2)
    pooledData(cond).matrix = [d1(cond).matrix; d2(cond).matrix];
    if size(d1,2) > 108
        pooledData(cond).gpfaTraj = [d1(cond).gpfaTraj; d2(cond).gpfaTraj];
    end
end

end