function output = pairPointPosVel(kin,closeFarther,distCriterion)
%% OVERVIEW

% This function pairs points in time by fulilling a set of criterion for
% distances in position and velocity space. closeFarther is a set of flags
% indicating whether, for position and velocity, the point should be
% further or closer together than the  distances in distCriterion.

% closeFarther ranges over three values. 1 indicates points should be
% closer than, -1 indicates points should be further than, and 0 indicates
% the distance does not matter. 

%% Pair point.

% Loop over condition, collect vars.
for cond = 1:size(kin,2)
    pos(cond).matrix = kin(cond).matrix(1:end-1,:);
    vel(cond).matrix = diff(kin(cond).matrix);
end

% Estabish a matrix for referencing indices. 
condRef = [];
timeRef =[];
for cond = 1:1:size(kin,2)
    condRef = [condRef repmat(cond,1,size(pos(cond).matrix,1))];
    timeRef = [timeRef 1:size(pos(cond).matrix,1)];
end
% condRef = kron(1:size(kin,2),ones(1,size(pos(1).matrix,1))).';
% timeRef = kron(ones(1,size(kin,2)),1:size(pos(1).matrix,1)).';

% Get the pairwise distance matrices.
pos2 = vertcat(pos.matrix);
vel2 = vertcat(vel.matrix);
posDists = squareform(pdist(pos2));
velDists = squareform(pdist(vel2));

% Loop over collected time points. 
for point = 1:size(posDists,1)
    % Assign the ID of that point.
    output(condRef(point),timeRef(point)).ID = [condRef(point) timeRef(point)];
    output(condRef(point),timeRef(point)).kin = ...
        [pos2(point,:) vel2(point,:)];
    % Find points that fulfill the position criterion.
    posVec = posDists(:,point);
    if closeFarther(1) == 1
        useIndsPos = find(posVec < distCriterion(1));
    elseif closeFarther(1) == -1
        useIndsPos = find(posVec > distCriterion(1));
    elseif closeFarther(1) == 0
        useIndsPos = 1:length(posVec);
    end
    % Find points that fulfill the velocity criterion.
    velVec = velDists(:,point);
    if closeFarther(2) == 1
        useIndsVel = find(velVec < distCriterion(2));
    elseif closeFarther(2) == -1
        useIndsVel = find(velVec > distCriterion(2));
    elseif closeFarther(2) == 0
        useIndsVel = 1:length(velVec);
    end
    % Combine the two, assign points.
    useInds = intersect(useIndsVel,useIndsPos);
    % Assign pairs.
    for cond = 1:size(kin,2)
        timeInds = useInds(find(condRef(useInds) == cond));
        output(condRef(point),timeRef(point)).pairPoints(cond).timepoints = ...
            timeRef(timeInds);
%         output(condRef(point),timeRef(point)).pairPoints(cond).kin = ...
%             mean([pos(cond).matrix(timeRef(timeInds),:) ...
%             vel(cond).matrix(timeRef(timeInds),:)],1);
    end
end






