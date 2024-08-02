function map = fitMap(predictor,predicted,lag)
%% OVERVIEW

% This function just automates fitting linear regression from one set of
% timeseries to another set of time series with lag. This function further
% detects and filters out nans from the fitting.

%% Fit map.

% Incorporate lag.
for cond = 1:size(predicted,2)
    predictor(cond).matrix = predictor(cond).matrix(:,1:end-lag);
    predicted(cond).matrix = predicted(cond).matrix(:,1+lag:end);
end

% Filter out nans.
for cond = 1:size(predicted,2)
    inds = unique([find(isnan(predictor(cond).matrix(1,:))) ...
        find(isnan(predicted(cond).matrix(1,:)))]);
    predictor(cond).matrix(:,inds) = [];
    predicted(cond).matrix(:,inds) = [];
end

% Fit the map.
map = [predicted().matrix]*pinv([predictor().matrix]);

end