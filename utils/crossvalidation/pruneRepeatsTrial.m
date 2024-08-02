function data = pruneRepeatsTrial(data)
%% OVERVIEW

% This function prunes out the trials that involve a visual distractor
% in the monkeys vision, but are essentially repeats of curved conditions.

%% Prune repeats.

% Get the conditions to assemble.
conds = unique([data().condNum]);

% Get the trials corresponding to distractor conds.
inds = ismember([data().condNum],find(mod(conds,3) == 0));

% Prune.
data(inds) = [];

end