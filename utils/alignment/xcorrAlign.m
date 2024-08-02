function output = xcorrAlign(data,filterSize,maxLag,ref)
%% OVERVIEW

% This function attempts to find an alignment between multiple instances of
% the "same" data by sliding to a reference condition. This is done by
% using finding a lag that maximizes the inner product of another condition
% with the reference condition. 

%% Align datasets.

% Filter datasets.
for cond = 1:size(data,2)
    for neuron = 1:size(data(cond).matrix,1)
        dataUse(cond).matrix(neuron,:) = imgaussfilt( ...
            data(cond).matrix(neuron,:),filterSize,'padding','symmetric', ...
            'filterSize',1+filterSize*10);
    end
end

% Loop over conditions, attempt to align datasets.
for cond = 1:size(data,2)
    if cond == ref
        output.aligned(cond).time = ...
            1:size(data(cond).matrix,2);
        output.aligned(cond).matrix = data(cond).matrix;
        ref
    else
        % Assign "time-points". 
        timePoints = 1:max([size(data(cond).matrix,2) size(data(ref).matrix,2)]);
        % Assign the correlations.
        output.corr(cond).vec = zeros(maxLag*2+1,1);
        % Loop over lags.
        for lag = -maxLag:maxLag
            if lag <= 0
                useRef = dataUse(ref).matrix(:,-lag+1:end);
                useAlign = dataUse(cond).matrix;
            elseif lag > 0
                useRef = dataUse(ref).matrix;
                useAlign = dataUse(cond).matrix(:,lag:end);
            end
            useRef = useRef(:,1:min([size(useRef,2) size(useAlign,2)]));
            useAlign = useAlign(:,1:min([size(useRef,2) size(useAlign,2)]));
            % Measure similarity.
            score = ...
                useRef(:).'*useAlign(:)/(norm(useRef(:),2)*norm(useAlign(:)));
            %score = nanmean(diag(corr(useAlign.',useRef.')));
            % Assign the inner product.
            output.corr(cond).vec(1+lag+maxLag) = score;
        end
        % Get the maximum alignment.
        [~,ind] = max(output.corr(cond).vec);
        % Assign time.
        output.aligned(cond).time = ...
            [1:size(data(cond).matrix,2)] ...
            -(ind-1-maxLag);
        output.aligned(cond).matrix = data(cond).matrix;
    end
end


end

