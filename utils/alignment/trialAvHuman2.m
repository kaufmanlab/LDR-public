function FRs = trialAvHuman2(trialStruct,smooth)
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
    % Out.
    count = count+1;
    % Find indices that correspond to that target.
    inds = find([trialStruct.kinematics(useInds).angle]==targs);
    inds2 = find([trialStruct.kinematics(useInds).out]==1);
    inds = intersect(inds,inds2);
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
                2,'padding','symmetric','filterSize',31);
        end
    end
    FRs.neuralActivity(count).matrix = neuralAct*50;
    FRs.kinematics(count).matrix = kin;
    FRs.kinematics(count).angle = trialStruct.kinematics(useInds(ind)).angle;
    FRs.kinematics(count).target = trialStruct.kinematics(useInds(ind)).target;
    FRs.kinematics(count).move = trialStruct.kinematics(useInds(ind)).move;
    FRs.kinematics(count).out = 1;
    % Back.
    count = count+1;
    % Find indices that correspond to that target.
    inds = find([trialStruct.kinematics(useInds).angle]==targs);
    inds2 = find([trialStruct.kinematics(useInds).out]==0);
    inds = intersect(inds+1,inds2);
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
                2,'padding','symmetric','filterSize',31);
        end
    end
    FRs.neuralActivity(count).matrix = neuralAct*50;
    FRs.kinematics(count).matrix = kin;
    FRs.kinematics(count).angle = FRs.kinematics(count-1).angle;
    FRs.kinematics(count).target = FRs.kinematics(count-1).target;
    FRs.kinematics(count).move = trialStruct.kinematics(useInds(ind)).move;
    FRs.kinematics(count).out = 0;
end

end