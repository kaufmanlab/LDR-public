function [trialKinematics,predictedKinematics] = ...
    predictKinUsingTBFandDataCondStructureSingletrial(kinematics,trialKinematics, ...
    leaveOut,data,trialData,L2,factorDim)
%% OVERVIEW

% This function attempts to predict kinematics from the corresponding
% neural activity, in particular using the loading matrix of TBF as 
% regressors. This attempts to linearly link the dynamics of motor cortex
% with produced reaches. This function uses linear regression to link the
% two. This function extends to single trials by fitting the map on trial
% averaged and then using it one single trials.

% The parameter L2 specifies an L2 penality on the regressors. 

%% Predict the held-out data.

% Prepare the kinematics.
predictedKinematics = trialKinematics;
kinFxns = kinFuncEval(-350:10:850,'position');
for cond = 1:size(kinematics,2)
    kinematics(cond).matrix = [kinematics(cond).X kinematics(cond).Y];
    kinematics(cond).matrix = ...
        pinv(kinFxns)*kinematics(cond).matrix;
    kinematics(cond).matrix = kinematics(cond).matrix(:);
end

% Create partitions of the data.
[train,test] = splitTrials([data().condNum],leaveOut);

% Get the trial condition numbers.
condInds = [trialData().condNum];
condNums = [data().condNum];

% Loop through the partition.
for partition = 1:size(train,1)
    % Get loading matrices.
    [basisFxns,loadings] = getBasisFxns(data(find(train(partition,:))),factorDim);
    for cond = 1:size(loadings,2)
        loadings(cond).matrix = loadings(cond).matrix(:);
    end
    heldInData = [loadings().matrix];
    heldInData = [heldInData; ones(1,size(heldInData,2))];
    % Get the held-in kinematics.
    heldInKin = [kinematics(find(train(partition,:))).matrix];
    % Get the map.
    linearMap = heldInKin*heldInData.'*pinv(heldInData*heldInData.' ...
        +L2*eye(size(heldInData,1)));
    % Test new trials.
    predConds = find(test(partition,:));
    predTrials = find(ismember(condInds,condNums(predConds)));
    for trial = 1:length(predTrials)
        % Get the loading matrix for the held-out condition.
        testLoading = 100*rebin(trialData(predTrials(trial)).matrix,10)*pinv(basisFxns.');
        testLoading = testLoading(:);
        % Get the prediction of the data.
        predictedKinematics(predTrials(trial)).matrix = ...
            linearMap(:,1:end-1)*testLoading+linearMap(:,end);
    end
end

% Decompress.
kinFxnsVel = kinFuncEval(-350:10:850,'velocity');
for cond = 1:size(trialKinematics,2)
    predictedKinematics(cond).velocity = kinFxnsVel*reshape( ...
        predictedKinematics(cond).matrix, ...
        [2 2]);
    predictedKinematics(cond).position = kinFxns*reshape( ...
        predictedKinematics(cond).matrix, ...
        [2 2]);
    % For analysis.
    predictedKinematics(cond).coeffs = predictedKinematics(cond).matrix;
    trialKinematics(cond).coeffs = ...
        pinv(kinFxns)*[trialKinematics(cond).X trialKinematics(cond).Y];
    trialKinematics(cond).coeffs = trialKinematics(cond).coeffs(:);
    % For quantification.
    trialKinematics(cond).matrix = [trialKinematics(cond).X trialKinematics(cond).Y];
end

end
