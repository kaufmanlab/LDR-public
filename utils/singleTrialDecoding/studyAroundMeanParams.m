function [correlations,pVals,coeffs,predCoeffs] = studyAroundMeanParams(kinematics,predictedKinematics)
%% OVERVIEW

% This function quantifies to what extent single-trial decoding acts as a
% classifier, and to what extent it acts as a continous decoder. It does by
% for each condition quantifying the correlations between the coefficients
% for the kinematic basis for the decoded and actual kinematics, thus
% subtracting off the average parameters and studying the variation around
% the condition average.

%% Quantify the correlations.

% Subtract off trial-averages.
kinForAv = kinematics;
predictedKinForAv = predictedKinematics;
for trial = 1:size(kinematics,2)
    predictedKinematics(trial).coeffs = ...
        predictedKinematics(trial).coeffs-mean( ...
        [predictedKinForAv([predictedKinematics().condNum] ...
        == predictedKinForAv(trial).condNum).coeffs],2);
    kinematics(trial).coeffs = ...
        kinematics(trial).coeffs-mean( ...
        [kinForAv([kinematics().condNum] ...
        == kinForAv(trial).condNum).coeffs],2);
end
[correlations,pVals] = corr( ...
        [kinematics().coeffs].', ...
        [predictedKinematics().coeffs].');
correlations = diag(correlations);
pVals = diag(pVals);
coeffs = [kinematics().coeffs].';
predCoeffs = [predictedKinematics().coeffs].';

end

