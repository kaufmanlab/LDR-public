function useEigs = bioEigs(useEigs)
%% OVERVIEW

% This function correct super-critical or stable eigenvalues to eigenvalues
% with a half life of ten seconds.

%% Correct eigs.

for i = 1:length(useEigs)
    if abs(useEigs(i)) > 0.9993
        useEigs(i) = 0.9993*useEigs(i)/abs(useEigs(i));
    end
end

end

