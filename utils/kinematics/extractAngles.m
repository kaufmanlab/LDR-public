function angs = extractAngles(Kinematics)
%% OVERVIEW

% This function takes kinematics and converts them to reach angles.

%% Extract angles.

angs = zeros(2,size(Kinematics,2));
for cond = 1:size(Kinematics,2)
    angs(:,cond) =  ...
        [(360/(2*pi))*atan2(Kinematics(cond).Y(end),Kinematics(cond).X(end)); ...
        cond];
end

end