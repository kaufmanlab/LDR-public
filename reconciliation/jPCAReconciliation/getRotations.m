function [rotations,reconstruction,varExplained,R2s] = getRotations(loadings,basisFxns,planeNum)
%% OVERVIEW

% This function takes in data, and uses jPCA to identify the rotations in
% the data. To allow reconstruction of the data, the mean subtraction and
% soft-normalization are turned off in jPCA and pre-applied. 

%% Identify the rotations in the data. 

% Reconstruct the data.
data = loadings;
for cond = 1:size(loadings,2)
    data(cond).matrix = loadings(cond).matrix*basisFxns.';
end

% Identify the mean response in the data.
meanResponse = 0;
for cond = 1:size(data,2)
    meanResponse = meanResponse + data(cond).matrix/size(data,2);
end

% Create an array for the usable data.
useData = data;
for cond = 1:size(data,2)
    useData(cond).A = data(cond).matrix.' - meanResponse.';
end

% Get the ranges for softnormalization.
softNorms = 1./(range(vertcat(useData().A))+5);
for cond = 1:size(data,2)
    useData(cond).A = softNorms.*useData(cond).A;
    useData(cond).A = useData(cond).A(20:40,:);
end

% Run jPCA and get rid of the annoying pop-ups.
[~,summary] = jPCA(useData);
R2s = [summary.R2_Mskew_2D summary.R2_Mbest_2D];
close all;

% Assign the rotations.
rotations = data;
for cond = 1:size(data,2)
    rotations(cond).matrix = ...
        summary.jPCs_highD(:,1:planeNum*2).'* ...
        (softNorms.'.*(data(cond).matrix-meanResponse));
    rotations(cond).jPCs = summary.jPCs_highD;
    rotations(cond).softNorms = softNorms;
end

% Assign the reconstruction.
reconstruction = data;
for cond = 1:size(data,2)
    reconstruction(cond).matrix = ...
        (summary.jPCs_highD(:,1:planeNum*2) ...
        *summary.jPCs_highD(:,1:planeNum*2).'* ...
        (softNorms.'.*(data(cond).matrix-meanResponse)))./softNorms.';
    reconstruction(cond).total = reconstruction(cond).matrix + meanResponse;
    reconstruction(cond).mean = meanResponse;
end
varExplained = getVarExplained(data,reconstruction,'ind');

end