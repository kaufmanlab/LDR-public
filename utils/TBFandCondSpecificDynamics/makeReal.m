function eigVecs = makeReal(eigVecs)
%% OVERVIEW

% This function accepts eigenvectors that may be complex, and returns the
% same eigenvectors but projected into the reals.

%% PARAMETERS

% Really annoyingly, small numerical errors means that real vectors can
% actually be imaginary with a imaginary norm on the order of 10e-17. This
% parameters just sets a bound below which imaginary norms are zero.
tol = 10^(-5);

%% MAKE THE EIGENVECTORS REAL

% For each eigenvectors, check if there is an imaginary component, if there
% is reassign as a real vector.
for col = 1:size(eigVecs,2)
    if norm(imag(eigVecs(:,col)),2) > tol
        realVec1 = eigVecs(:,col) + eigVecs(:,col+1);
        realVec2 = complex(0,1)*(eigVecs(:,col) - eigVecs(:,col+1));
        eigVecs(:,col) = real(realVec1/norm(realVec1,2));
        eigVecs(:,col+1) = real(realVec2/norm(realVec2,2));
    end
end

end

