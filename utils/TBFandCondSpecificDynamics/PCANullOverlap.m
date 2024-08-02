function nullIndices = PCANullOverlap(data,subspaces,method)
%% OVERVIEW

% This function takes in data and uses it to generate the singular vectors
% of the data in state space. This is then used to generate a null
% distribution of the expected alignment index between conditions. The
% "data" for the null distribution projected into the subspace identified
% by subspaces, and then re-oriented into a randomly-drawn subspace. The
% alignment indices are then quantified.

%% Generate null.

% Generate the SVs.
[U,S,~] = svd([data().matrix],'econ');
space = U*S;

% Generate the null.
nullIndices = zeros(size(data,2));
for cond1 = 1:size(data,2)
    for cond2 = 1:size(data,2)
        if strcmp(method,'loading')
            newSpace1 = space*normrnd(0,1, ...
                size(subspaces(cond1).matrix,1), ...
                size(subspaces(cond1).matrix,2));
            newSpace1 = newSpace1./sum(newSpace1.^2).^0.5;
            if cond1 ~= cond2
                newSpace2 = space*normrnd(0,1, ...
                    size(subspaces(cond2).matrix,1), ...
                    size(subspaces(cond2).matrix,2));
                newSpace2 = newSpace2./sum(newSpace2.^2).^0.5;
            else
                newSpace2 = newSpace1;
            end
            nullIndices(cond1,cond2) = getAlignmentIndex(...
                orth(newSpace1)*pinv(subspaces(cond1).eigVecs)*data(cond1).matrix, ...
                newSpace2);
        elseif strcmp(method,'randDraw')
            newCond1 = space*normrnd(0,1,size(data(cond1).matrix,1), ...
                size(data(cond1).matrix,2));
            if cond1 ~= cond2
                newCond2 = space*normrnd(0,1,size(data(cond2).matrix,1), ...
                    size(data(cond2).matrix,2));
            else
                newCond2 = newCond1;
            end
            [compareSpace,~,~] = svd(newCond2,'econ');
            compareSpace = compareSpace(:,1:size(subspaces(cond2).matrix,2));
            nullIndices(cond1,cond2) = getAlignmentIndex(...
                newCond1,compareSpace);
        end
    end
end

end