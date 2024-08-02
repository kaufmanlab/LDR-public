function spikesNewBins = rebin(spikes,ratio)
%% OVERVIEW

% This function converts spike counts in one length time bin to another
% length time bin. To avoid workin in units, the new time bins are only
% given as a length ratio to the original time bins. To avoid frequency
% artifacts when converting from bins to spike times and to account for the
% lack of knowledge as to "where" the spike was in the original time bin,
% spikes are converted to the earliest edge of their original time bin then
% peturbed by a peturbation uniformly distributed across the time bin,
% jittering the spike. 

%% Rebin spikes.

% Convert from sparse and logical just in case.
spikes = full(double(spikes));

% Get the bins.
newBins = 0:ceil(size(spikes,2)/ratio);
newBins = newBins*ratio;

% Generate the new bins.
spikesNewBins = zeros(size(spikes,1),length(newBins)-1);

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
    % Assign the new spikes.
    [spikesNewBins(neuron,:),~] = histcounts(spikeTime,newBins);
end

end

