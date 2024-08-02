function data = pruneRepeats(data)
%% OVERVIEW

% This function prunes out the conditions that involve a visual distractor
% in the monkeys vision, but are essentially repeats of curved conditions.

%% Prune repeats.


data(3:3:end) = [];

end