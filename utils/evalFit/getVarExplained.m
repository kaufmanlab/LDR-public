function VE = getVarExplained(data,approx,flag)
%% OVERVIEW

% This function takes data and approximations to the data and returns the 
% quality of the fit as quantified by the variance explained. If the flag
% is 'pooled', then a scalar is returned. If the flag is 'ind', then an
% array of variance explained for single conditions is returned.

%% Quantify the variance explained.

% Remove parts of the data not approximated, assuming that approximation is
% uniform across neurons but not time.
for cond = 1:size(data,2)
    inds = find(isnan(approx(cond).matrix(1,:)));
    approx(cond).matrix(:,inds) = [];
    data(cond).matrix(:,inds) = [];
end

% Quantify the variance explained.
if strcmp(flag,'pooled')
    VE = 1 -sum(var([data().matrix] - [approx().matrix],0,2))/ ...
        sum(var([data().matrix],0,2));
elseif strcmp(flag,'ind')
    VE.array = zeros(1,size(data,2));
    for cond = 1:size(data,2)
        VE.array(cond) = 1 -sum(var(data(cond).matrix - approx(cond).matrix,0,2))/ ...
            sum(var(data(cond).matrix,0,2));
    end
    VE.array(VE.array < 0) = 0;
    VE.mean = mean(VE.array);
    VE.std = std(VE.array);
end

end

