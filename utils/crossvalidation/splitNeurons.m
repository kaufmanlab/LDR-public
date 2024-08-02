function [predictors,predicted] = splitNeurons(data,ratio)
%% OVERVIEW

% This function takes in population recordings, partitions the neurons in
% the dataset into predictor neurons and predicted neurons, and then
% returns two seperate datasets where the neurons have been partitioned.

%% Partition neurons.

% Get a permutation of neurons.
neuronInds = 1:size(data(1).matrix,1);
neuronInds = neuronInds(randperm(length(neuronInds)));

% Partition the neurons.
trainInds = neuronInds(1:round(length(neuronInds)*ratio));
testInds = neuronInds(1+round(length(neuronInds)*ratio):end);

% Assign arrays then partition.
predictors = data;
predicted = data;
for cond = 1:size(data,2)
    predictors(cond).matrix = data(cond).matrix(trainInds,:);
    predicted(cond).matrix = data(cond).matrix(testInds,:);
end

end