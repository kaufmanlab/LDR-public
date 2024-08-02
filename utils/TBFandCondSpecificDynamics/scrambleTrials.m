function trials = scrambleTrials(trials)
%% OVERVIEW

% This function scrambles the temporal indices of trials randomly via
% permutations, destroying the temporal characteristics of the dataset
% while maintaining the total spike counts and the spike correlations
% between neurons. 

%% Scramble trials. 

for trial = 1:size(trials,2)
    trials(trial).matrix = full(double(trials(trial).matrix));
    trials(trial).matrix = trials(trial).matrix(:,randperm(size(trials(trial).matrix,2)));
end

end