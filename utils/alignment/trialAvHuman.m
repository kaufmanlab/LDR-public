function FRs = trialAvHuman(trialStruct,smooth)
%% OVERVIEW

% This function takes in a data structure produced by Sliman et al, in
% particular of a participant viewing a center-out reaching task. This is
% divided up into trials that are aligned to movement onset. 

%% Get trials.


% Get the unique targets.
%useInds = find([trialStruct.kinematics.observed]);
useInds = 1:size(trialStruct.kinematics,2);
uniqueTargs = unique([trialStruct.kinematics(useInds).angle]);
% Loop over targets.
count = 0;
for targs = uniqueTargs
    count = count+1;
    % Find indices that correspond to that target.
    inds = find([trialStruct.kinematics(useInds).angle]==targs);
    % Loop over the indices.
    neuralAct = 0;
    kin = 0;
    targ = 0;
    for ind = inds
        neuralAct = neuralAct + trialStruct.neuralActivity(useInds(ind)).matrix/length(inds);
        kin = kin + trialStruct.kinematics(useInds(ind)).matrix/length(inds);
        targ = targ + trialStruct.kinematics(useInds(ind)).targetOverTime/length(inds);
    end
    targ(end-19:end,:) = repmat(targ(end-20,:),20,1);
    % Smooth if indicated.
    if smooth
        for neuron = 1:size(neuralAct,1)
            neuralAct(neuron,:) = ...
                imgaussfilt(neuralAct(neuron,:), ...
                2,'padding','symmetric','filterSize',31);
        end
        for dim = 1:size(kin,2)
            kin(:,dim) = ...
                imgaussfilt(kin(:,dim), ...
                2,'padding','symmetric','filterSize',31);
        end
    end
    FRs.neuralActivity(count).matrix = neuralAct*50;
    FRs.kinematics(count).matrix = kin;
    FRs.kinematics(count).targOverTime = targ;
    FRs.kinematics(count).angle = trialStruct.kinematics(useInds(ind)).angle;
    FRs.kinematics(count).target = trialStruct.kinematics(useInds(ind)).target;
    FRs.kinematics(count).move = trialStruct.kinematics(useInds(ind)).move;
    FRs.kinematics(count).out = 1;
end

end