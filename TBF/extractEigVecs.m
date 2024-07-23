function eigVecs = extractEigVecs(PA,frequencies)
%% OVERVIEW

% This function takes the principal axes of rotational planes for each
% condition after alignment between conditions and extracts an eigen-basis
% that describes an equivalent ellipsoid in state space. The intution for
% this procedure is that while there are unbounded ways of describing an
% ellipsoid or rotational dynamics, the "proper eigenbasis" should be
% unit-norm in both vectors and the projection of activity should always be
% circular, so that any oblongness in state space arises from correlations
% between eigenvectors. This basis can be found by rotating the principal
% axes by 45 degrees. 

%% Extract the eigenbasis.

% Get the IDs of pairs.
ID = pairID(frequencies);

% Go through the straight conditions and align principal axes.
eigVecs = PA;
ang = 2*pi*45/360;
for cond = 1:size(PA,2)
    for el = unique(ID)
        inds = find(ID == el);
        if length(inds) == 2
            eigVecs(cond).matrix(:,inds) = ...
                PA(cond).matrix(:,inds) ...
                *[cos(ang) sin(ang); -sin(ang) cos(ang)];
        end
        eigVecs(cond).matrix(:,inds) = ...
            eigVecs(cond).matrix(:,inds) ...
            ./(sum(eigVecs(cond).matrix(:,inds).^2).^0.5);
    end
end

end