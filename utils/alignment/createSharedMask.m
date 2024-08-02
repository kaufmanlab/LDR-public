function mask = createSharedMask(masks)
%% OVERVIEW

% This function creates a mask that is shared among multiple datasets.

%% Create shared mask.

mask = logical(min(masks,[],2));

end