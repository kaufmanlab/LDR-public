function prediction = applyMap(predictor,map,lag)
%% OVERVIEW

% This function takes in predictors for times series, the map from the
% predictor to the predicted time series, and the lag between the two, 
% applies the transformation and returns the predictions.

%% Predict.

% Apply the transformation.
prediction = predictor;
for cond = 1:size(predictor,2)
    prediction(cond).matrix = map*predictor(cond).matrix(:,1:end-lag);
end

end