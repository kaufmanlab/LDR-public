function [train,test] = splitTrials(trials,ratio)
%% OVERVIEW

% This function splits trials indices into train-test partitions based on
% the specified ratio. If the ratio is not a divisor of the number of
% trials, the function rejects everything. If ratio is
% "leaveOneOut", then the ratio is set to the number of trials.

%% Partition the trials. 

% HACK FOR TRAINING ON LESS THAN DATA THAN TESTING.
flip = false;
if ratio < 1
    flip = true;
    ratio = 1/ratio;
end

% Reject if the ratio is not a divisor.
if ~strcmp(ratio,'leaveOneOut')
    if mod(length(unique(trials))/ratio,1) > 0
        disp('Specify integer divisor');
        return
    end
end

% Modify to leave one out if requested.
if strcmp(ratio,'leaveOneOut')
    ratio = length(trials);
end
    
% Assign the arrays.
train = zeros(ratio,length(trials));
test = zeros(ratio,length(trials));

% Randomly permute the unique indices.
uniqueTrials = unique(trials);
inds = 1:length(uniqueTrials);
newInds = randperm(length(uniqueTrials));

% Assign the partitions.
for split = 1:ratio
    testInds = inds((1:length(uniqueTrials)/ratio)+(split-1)*length(uniqueTrials)/ratio);
    trainInds = inds(~ismember(inds,testInds));
    train(split,:) = ismember(trials,uniqueTrials(newInds(trainInds)));
    test(split,:) = ismember(trials,uniqueTrials(newInds(testInds)));
end
train = logical(train);
test = logical(test);

% HACK FOR TRAINING ON LESS THAN DATA THAN TESTING.
if flip
    copyTrain = train;
    copyTest = test;
    test = copyTrain;
    train = copyTest;
end

% HACK FOR COMPARING TO RNN.
if ratio == 1
    train = train*0 + 1;
    test = test*0 + 1;
end

end

