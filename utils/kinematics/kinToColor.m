function color = kinToColor(X,Y)
%% OVERVIEW

% This function takes in the end-points of a reach, and returns an RGB
% vector that dictates the color of the corresponding PSTH. 

%% CREATE THE COLOR

% Get the angle.
angle = atan2(X,Y);

% Assign the color.
color = [(1+sin(angle+pi/2))/2 ...
    (1+sin(angle))/2 ...
    (1+sin(-angle-pi/2))/2]*0.85;

end

