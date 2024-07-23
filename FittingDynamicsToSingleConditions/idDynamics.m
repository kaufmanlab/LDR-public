function [subspaces,eigVals,params,approx,sim,basisFxns] = idDynamics(data,dimUse)
%% OVERVIEW

% This function takes in data locked to movement onset, along with a
% parameter expressing the dimensionality to use, and returns for each
% condition independently the eigenvectors in which the dynamics reside,
% along with the eigenvalues of those dynamics, the parameters described by
% those eigenvalues, and the approximation to the data provided. 

%% Find the dynamics describing each conditions.

% For each condition, find dynamics. 
subspaces = data;
eigVals = zeros(dimUse,size(data,2));
approx = data;
basisFxns = zeros(size(data(1).matrix,2),dimUse,size(data,2));
params.param = zeros(dimUse,3,size(data,2));
params.functionalForm = zeros(dimUse,3,size(data,2));
params.meaning = ...
    ["Frequency (cycles/sec)" "Phase (degrees)" "Half-life (sec)"];
for cond = 1:size(data,2)
    % Fit dynamics to the data in the low-dimensional space. 
    dynMatrix = (data(cond).matrix(:,21:end)) ...
        *pinv(data(cond).matrix(:,20:end-1));
    % Solve RRR.
    [subspace,~,~] = svd(dynMatrix*data(cond).matrix(:,20:end),'econ');
    subspace = subspace(:,1:dimUse);
    dynMatrix = (subspace*subspace.')*dynMatrix;
    % Eigendecompose the dynamics.
    [E,e] = eigs(dynMatrix,dimUse);
    e = bioEigs(diag(e));
    [~,I] = sort(atan2(abs(imag(e)),real(e))/(2*pi));
    e = e(I);
    E = E(:,I);
    eigVals(:,cond) = e;
    % Find the pairs of eigenvalues.
    ID = pairID(eigVals(:,cond));
    % Assign the parameters.
    for el = 1:length(unique(ID))
        inds = find(ID == el);
        useEig = e(inds(1));
        [f,t,frequency,halflife] = extractParams(useEig);
        params.functionalForm(inds(1),:,cond) = [f 0 t];
        params.param(inds(1),:,cond) = [frequency 0 halflife];
        if length(inds) == 2
            params.functionalForm(inds(2),:,cond) = [f -pi/2 t];
            params.param(inds(2),:,cond) = [frequency -90 halflife];
        end
    end
    % Assign the basis functions.
    timeStamp = [nan(19,1); (0:101).'];
    for col = 1:size(basisFxns,2)
        basisFxns(:,col,cond) = ...
            exp(-params.functionalForm(col,3,cond)*timeStamp) ...
            .*cos(2*pi*params.functionalForm(col,1,cond)*timeStamp ...
            + params.functionalForm(col,2,cond));
    end
    % Assign the approximation.
    approx(cond).matrix = data(cond).matrix;
    sim(cond).matrix = data(cond).matrix;
    for t = 21:121
        sim(cond).matrix(:,t) = ...
            dynMatrix*sim(cond).matrix(:,t-1);
    end
    sim(cond).matrix(:,1:19) = sim(cond).matrix(:,1:19)*nan;
    % Assign the subspace.
    subspaces(cond).matrix = E;
end
    
end

