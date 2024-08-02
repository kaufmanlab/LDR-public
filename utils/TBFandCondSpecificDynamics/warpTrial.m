function spikes = warpTrial(spikes,warpFrom,warpTo)
%% OVERVIEW

% This function converts spike counts to raw time by assigning the spike
% time to the front edge of each 1 ms bin then jittering the spike with a
% uniform distribution over the bin. This is done so that the spike times
% may then be linearly warped such that all reaches last the same length of
% time. Spikes are then re-binned in the same units.

%% Rebin spikes.

% Assign the lag.
lag = 100;

% Convert from sparse and logical just in case.
spikes = full(double(spikes));

% For each neuron, peturb then re-bin.
for neuron = 1:size(spikes,1)
    % Initialize the count.
    spikeCounter = 0;
    % Get the total number of spikes.
    spikeNum = sum(spikes(neuron,:));
    % Assign an array for spike times.
    spikeTime = zeros(1,spikeNum);
    % Get the indices that have a spike in them. 
    spikeInds = find(spikes(neuron,:) > 0);
    % Run through the spikes.
    for spikeInd = spikeInds
        % Get the number of spikes.
        emitted = spikes(neuron,spikeInd);
        % Run through the individual spikes, append to array.
        for spike = 1:emitted
            spikeCounter = spikeCounter + 1;
            spikeTime(spikeCounter) = spikeInd-1+rand;
        end
    end
    % Put the spikes in the correct time frame, with additional lag built
    % for correction.
    spikeTime = spikeTime-350+lag;
    % Add a jitter to the reach end.
    jitter = 0;
    % Find spikes that happen during reaching and after the reach.
    reachSpikeInd = find(spikeTime > 0 & spikeTime < warpFrom+jitter);
    postReachSpikeInd = find(spikeTime >= warpFrom+jitter);
    % Warp the spikes.
    spikeTime(reachSpikeInd) = ...
        (warpTo)*spikeTime(reachSpikeInd)/(warpFrom+jitter);
    %spikeTime(postReachSpikeInd) = (850+lag-(warpTo+jitter))*...
    %    (spikeTime(postReachSpikeInd) - warpFrom+jitter) ...
    %    /(850+lag-(warpFrom+jitter))+(warpTo+jitter);
    % Assign the new spikes.
    spikeTime(postReachSpikeInd) = spikeTime(postReachSpikeInd) ...
        -(warpFrom+jitter)+(warpTo+jitter);
    spikeTime = spikeTime+350-lag;
    [spikes(neuron,:),~] = histcounts(spikeTime,0:1201);
end

end

