function [train,test] = splitTrialsAngles(trials,ratio,Kinematics)
%% OVERVIEW

% This function splits trials indices into train-test partitions based on
% the specified ratio. If the ratio is not a divisor of the number of
% trials, the function rejects everything. If ratio is
% "leaveOneOut", then the ratio is set to the number of trials.

% Unlike a previous version, this version first calculates the angles of
% reaches and sorts them, then leaves out fractions of the reach space.

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

% Assign angles. 
angles = zeros(1,size(Kinematics,2));
for cond = 1:length(angles)
    if isfield(Kinematics(1),'Y')
        angles(cond) = atan2(Kinematics(cond).Y(end),Kinematics(cond).X(end));
    else
        angles(cond) = Kinematics(angle);
    end
end
[~,newInds] = sort(angles);

% Assign the partitions.
for split = 1:ratio
    testInds = newInds((1:length(newInds)/ratio)+(split-1)*length(newInds)/ratio);
    trainInds = newInds(~ismember(newInds,testInds));
    train(split,trainInds) = 1;% = ismember(trials,trainInds);
    test(split,testInds) = 1;%ismember(trials,testInds);
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

end


