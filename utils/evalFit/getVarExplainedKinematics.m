function VE = getVarExplainedKinematics(kinematics,approx)
%% OVERVIEW

% This function takes kinematics and approximations to the kinematics and  
% return the quality of the fit as quantified by the variance explained for
% the position and the velocity.

%% Quantify the variance explained.

% Quantify the variance explained.
VE.pos.array = zeros(1,size(kinematics,2));
for cond = 1:size(kinematics,2)
    VE.pos.array(cond) = 1 -sum(var(kinematics(cond).matrix - approx(cond).position,0,1))/ ...
        sum(var(kinematics(cond).matrix,0,1));
end
VE.pos.array(VE.pos.array < 0) = 0;
VE.pos.mean = mean(VE.pos.array);
VE.pos.std = std(VE.pos.array);

VE.vel.array = zeros(1,size(kinematics,2));
for cond = 1:size(kinematics,2)
    VE.vel.array(cond) = 1 -sum(var(diff(kinematics(cond).matrix,1,1) - approx(cond).velocity(2:end,:),0,1))/ ...
        sum(var(diff(kinematics(cond).matrix,1,1),0,1));
end
VE.vel.array(VE.vel.array < 0) = 0;
VE.vel.mean = mean(VE.vel.array);
VE.vel.std = std(VE.vel.array);
    
end

