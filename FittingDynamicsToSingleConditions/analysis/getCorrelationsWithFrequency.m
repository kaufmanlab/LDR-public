function output = getCorrelationsWithFrequency
%% OVERVIEW

% This function checks whether there exist correlations in the frequency
% and half life of single-condition dynamics with kinematics. In
% particular, we check the reach duration, curvature, direction.

%% Compare trials to the collected distribution.

% Load the data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);

% Load the results.
load('condSpecificDynamics');
load('warpedData');

% For each monkey analyze.
for monkey = 1:size(ShenoyMonkeyData,2)
    output(monkey).M1 = analysis(condSpecificDynamics(monkey).M1.eigVals, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        [warpedData(monkey).M1.oldDuration]);
    output(monkey).PMd = analysis(condSpecificDynamics(monkey).PMd.eigVals, ...
        ShenoyMonkeyData(monkey).Kinematics, ...
        [warpedData(monkey).PMd.oldDuration]);
end

end

%% Subfunction for quantifying correlations. 

function output = analysis(eigVals,kinematics,durations)
kinematics = pruneRepeats(kinematics);
% Assemble various representations of the eigenvalues.
realImagRep = [real(eigVals.') imag(eigVals.')];
%angleMagRep = [abs(eigVals.') atan2(imag(eigVals.'),real(eigVals.'))];
output.allRepEigVals = realImagRep;
output.allRepEigVals = [output.allRepEigVals ones(size(output.allRepEigVals,1),1)];
output.allRepEigValsLabels = ...
    [repmat("real",1,size(eigVals,1)) repmat("imag",1,size(eigVals,1))];
% Assemble kinematic functions.
kinFuncs = kinFuncEval(-350:10:850,'position');
% Assemble kinematic correlates.
kinCorrs = zeros(size(kinematics,2),5);
for cond = 1:size(kinematics,2)
    params = pinv(kinFuncs)*[kinematics(cond).X kinematics(cond).Y];
    kinCorrs(cond,:) = [params(1,:) params(2,:) durations(cond)];
end
output.kinParams = kinCorrs;
output.kinParams = [output.kinParams ones(size(output.kinParams,1),1)];
output.kinLabels = ["x-extent" "y-extent" "x-curve" "y-curve" "duration"];
% Quantify correlations.
[output.corrs,output.pVals] = corr(output.kinParams,output.allRepEigVals);
% Quantify predictability.
output.allRepEigValsHat = output.allRepEigVals;
output.kinParamsHat = output.kinParams;
holdSplit = 72;
[train,test] = splitTrials(1:72,holdSplit);
for fold = 1:holdSplit
    % For matrix from neural to kinematics.
    neuralToKin = pinv(output.allRepEigVals(find(train(fold,:)),:)) ...
        *output.kinParams(find(train(fold,:)),:);
    % For matrix from kinematics to neuural. 
    kinToNeural = pinv(output.kinParams(find(train(fold,:)),:)) ...
        *output.allRepEigVals(find(train(fold,:)),:);
    % Go over training conditions.
    for cond = find(test(fold,:))
        output.allRepEigValsHat(cond,:) = output.kinParams(cond,:)*kinToNeural;
        output.kinParamsHat(cond,:) = output.allRepEigVals(cond,:)*neuralToKin;
    end
end
% Calculate R2.
output.allRepEigValsVarExpl = ...
    1-var(output.allRepEigValsHat - output.allRepEigVals,0,1) ...
    ./var(output.allRepEigVals,0,1);
output.allRepEigValsVarExpl(output.allRepEigValsVarExpl < 0) = 0;
output.kinParamsExpl = ...
    1-var(output.kinParamsHat - output.kinParams,0,1) ...
    ./var(output.kinParams,0,1);
output.kinParamsExpl(output.kinParamsExpl < 0) = 0;
end

