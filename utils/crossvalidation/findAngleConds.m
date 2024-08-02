function inds = findAngleConds(kinematics,angles)
%% OVERVIEW

% This function finds conditions that most closely align to a set of
% desired angles. 

%% Find minimizers.

% Assign indices.
inds = angles*0;

% Assign the end points.
endPoint = zeros(size(kinematics,2),2);
for cond = 1:size(kinematics,2)
    endPoint(cond,:) = [kinematics(cond).matrix(end,2) ...
        kinematics(cond).matrix(end,1)];
    endPoint(cond,:) = endPoint(cond,:)/norm(endPoint(cond,:),2);
end

% Assign the ideal end points.
angles = repmat(angles.',1,2);
for ang = 1:size(angles,1)
    angles(ang,:) = [cos(angles(ang,1)) sin(angles(ang,2))];
end

% Assign the nearest.
sims = endPoint*angles.';
for ang = 1:size(angles,1)
    [~,ind] = max(sims(:,ang));
    inds(ang) = ind;
end

end