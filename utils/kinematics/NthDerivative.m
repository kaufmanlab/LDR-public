function output = NthDerivative(output,n,type )
%% OVERVIEW

% This function differentiates or integrates a set of signals in output
% based on the provided n, where n > 0 is the nth derivative, n = 0 is the
% signal, and n < 0 is the nth integral (antiderivative) of the signal.

%% Diff or integrate.

if abs(n) > 2
    disp("Don't")
    return
end

if strcmp(type,'raw')
    if n > 0
        for cond = 1:size(output,2)
            output(cond).X = diff(output(cond).X,n,1);
            output(cond).X = [zeros(n,1); ...
                output(cond).X];
            output(cond).Y = diff(output(cond).Y,n,1);
            output(cond).Y = [zeros(n,1); ...
                output(cond).Y];
        end
    elseif n < 0
        for cond = 1:size(output,2)
            for op = 1:abs(n)
                output(cond).X = cumsum(output(cond).X);
                output(cond).Y = cumsum(output(cond).Y);
            end
        end
    end
else
    if n > 0
        for cond = 1:size(output,2)
            output(cond).matrix = diff(output(cond).matrix,n,1);
            output(cond).X = [zeros(n,size(output(cond).matrix,2)); ...
                output(cond).matrix];
        end
    elseif n < 0
        for cond = 1:size(output,2)
            for op = 1:abs(n)
                output(cond).matrix = cumsum(output(cond).matrix);
            end
        end
    end
end

end