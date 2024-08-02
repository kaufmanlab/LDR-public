function trials = preprocessTrials(trials)
%% OVERVIEW

% This function pre-processes trials binned at 1 ms resolution by first
% smoothing the trials with a 20 ms Gaussian kernel then down-sampling to
% 10 ms resolution.

%% Pre-process trials.

% For each trial, smooth then down-sample.
for trial = 1:size(trials,2)
    trialArray = full(double(trials(trial).matrix));
    for neuron = 1:size(trialArray)
        trialArray(neuron,:) = imgaussfilt(1000*trialArray(neuron,:), ...
            20,'FilterSize',111,'Padding','Symmetric');
    end
    trials(trial).matrix = trialArray(:,1:10:end);
end

end

