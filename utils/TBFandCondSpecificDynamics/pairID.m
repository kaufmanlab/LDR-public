function ID = pairID(eigVals)
%% OVERVIEW

% This function just identifies which pairs of eigenvalues are complex
% conjugates.

%% ID pairs.

resps = complex(real(eigVals),abs(imag(eigVals)));
[~,~,ID] = unique(resps);
if size(ID,2) == 1 
    ID = ID.';
end

end

