function [projection,dimension,AUC,pVal,varExplained] = ...
    identifyCIS(data,loadings,basisFxns,params,svNum)
%% OVERVIEW

% This function identifes the dimension of maximum seperation between the
% fixed points and the initial points, which is hypothesized to underlying
% the condition-invariant signal. This is done using Fisher Discriminant
% Analysis to find the dimension of maximum discriminability.

%% Identify CIS.

% Find the singular vectors of the dataset.
[regularize,~,~] = svd([data().matrix],'econ');
regularize = regularize(:,1:svNum);

% Get the fixed points.
ind = find(params(:,1) == 0);
temporalInds = find(~isnan(basisFxns(:,1)));
fpArray = zeros(size(data(1).matrix,1),size(data,2));
for cond = 1:size(data,2)
    fpArray(:,cond) = mean(loadings(cond).matrix(:,ind)*basisFxns(temporalInds,ind).',2);
end

% Get the initial states.
isArray = zeros(size(data(1).matrix,1),size(data,2));
for cond = 1:size(data,2)
    isArray(:,cond) = mean(data(cond).matrix(:,1:10),2);
end

% Regularize the array.
isArray = regularize.'*isArray;
fpArray = regularize.'*fpArray;

% Get the means.
fpMean = mean(fpArray,2);
isMean = mean(isArray,2);

% Get the covariances.
fpCov = cov(fpArray.');
isCov = cov(isArray.');

% Find the maximum discriminable dimension.
dimension = regularize*pinv(fpCov+isCov)*(fpMean-isMean);
dimension = dimension/norm(dimension,2);

% Get the projection.
projection = data;
for cond = 1:size(data,2)
    projection(cond).projected = dimension.'*projection(cond).matrix;
    projection(cond).matrix = dimension*projection(cond).projected;
end

% Get the ROC-AUC.
labels = [ones(1,size(data,2)) zeros(1,size(data,2))];
points = dimension.'*regularize*[isArray fpArray];
[~,~,~,AUC] = perfcurve(labels,points,0);
AUC = abs(AUC-0.5)+0.5;

% Control for significance.
drawNum = 1000;
control = zeros(1,drawNum);
for draw = 1:drawNum
    vec = normrnd(0,1,svNum,1);
    vec = vec/norm(vec,2);
    labels = [ones(1,size(data,2)) zeros(1,size(data,2))];
    points = vec.'*[isArray fpArray];
    [~,~,~,testauc] = perfcurve(labels,points,0);
    testauc = abs(testauc-0.5)+0.5;
    control(draw) = testauc;
end
pVal = length(find(control>AUC))/drawNum;

% Get the variance explained.
varExplained = getVarExplained(data,projection,'ind');

end

