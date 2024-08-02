function FRs = trialAvHumanTargetPursuit(trialStruct,holdOut,smooth)
%% OVERVIEW

% This function takes in a data structure produced by Sliman et al, in
% particular of a participant viewing a center-out reaching task. This is
% divided up into trials that are aligned to movement onset. 

%% Get trials.


% Get the unique targets.
%useInds = find([trialStruct.kinematics.observed]);
useInds = 1:size(trialStruct.kinematics,2);
% Assemble the to-from pairs.
pairs = [[trialStruct.kinematics.target]; [trialStruct.kinematics.from]];
% Use our pre-defined knowledge to know there are 20 pairings.
moveIndices = kmeans(pairs.',20).';
uniqueMoves = unique(moveIndices);
% Loop over targets.
count = 0;
for targs = uniqueMoves
    count = count+1;
    % Find indices that correspond to that target.
    inds = find(moveIndices==targs);
    % If any trials are indicated to be held out, permute then hold
    % out.
    if holdOut > 0
        inds = inds(randperm(length(inds)));
        holdOff = inds(end-holdOut+1:end);
        inds = inds(1:end-holdOut);
    end
    % Loop over the indices.
    neuralAct = 0;
    kin = 0;
    for ind = inds
        neuralAct = neuralAct + trialStruct.neuralActivity(useInds(ind)).matrix/length(inds);
        kin = kin + trialStruct.kinematics(useInds(ind)).matrix/length(inds);
    end
    % Smooth if indicated.
    if smooth
        for neuron = 1:size(neuralAct,1)
            neuralAct(neuron,:) = ...
                imgaussfilt(neuralAct(neuron,:), ...
                2.5,'padding','symmetric','filterSize',31);
        end
    end
    FRs.neuralActivity(count).matrix = neuralAct*50;
    FRs.kinematics(count).matrix = kin;
    FRs.kinematics(count).angle = trialStruct.kinematics(useInds(ind)).angle;
    FRs.kinematics(count).target = trialStruct.kinematics(useInds(ind)).target;
    FRs.kinematics(count).from = trialStruct.kinematics(useInds(ind)).from;
    FRs.kinematics(count).moveIndices = moveIndices(ind);
    if holdOut > 0
        for index = 1:length(holdOff)
            FRs.neuralActivity(count).heldOut(index).matrix = ...
                trialStruct.neuralActivity(useInds(holdOff(index))).matrix;
            FRs.kinematics(count).heldOut(index).matrix = ...
                trialStruct.kinematics(useInds(holdOff(index))).matrix;
        end
        for index = 1:length(inds)
            FRs.neuralActivity(count).heldIn(index).matrix = ...
                trialStruct.neuralActivity(useInds(inds(index))).matrix;
            FRs.kinematics(count).heldIn(index).matrix = ...
                trialStruct.kinematics(useInds(inds(index))).matrix;
        end
    end
end

end