function correction = produceCorrection(initialState,ID)
%% OVERVIEW

% This function take basis functions that, while is in an eigenbasis, are
% not aligned such that the initial state produces a sine/cosine pairs of
% unit magnitude, and produces a rotation and scaling matrix that corrects
% for this. 

% Just for clarity, this function technicallly produces the inverse of the
% transform to be applied to the loading matrix, assuming that the basis
% functions are being worked with implicitly. 

%% Correct initial state. 

% For each pair estimate the rotation needed.
correction = zeros(length(initialState),length(initialState));
for el = 1:length(unique(ID))
    inds = find(ID == el);
    state = initialState(inds);
    stateNorm = norm(state,2);
    if length(inds) == 2
        ang = atan2(state(2),state(1));
        rotMatrix = [cos(ang) -sin(ang); sin(ang) cos(ang)];
    else
        rotMatrix = sign(state);
    end
    correction(inds,inds) = stateNorm*rotMatrix;
end

end

