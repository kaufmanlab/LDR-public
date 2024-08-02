function lagged = createLags(mat,lags)
%% OVERVIEW

% This function takes in a matrix and produces multiple lagged versions of
% the matrix.

%% Lag

% Assign the output.
lagged = zeros(size(mat,1)-lags+1,size(mat,2)*lags);

% Populate.
for lag = 1:lags
    lagged(:,(1:size(mat,2))+(lag-1)*size(mat,2)) = mat(lag:end-lags+lag,:);
end

end