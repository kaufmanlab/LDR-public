function eigValsAv = averageEigVals(eigVals,flag)
%% OVERVIEW

% This function averages over eigenvalues. This has to be done in an
% intelligent manner (not just mean(eigVals)) because eigenvalues can
% either come in pairs or be independent on the real numbers, so just
% averaging over pairs where some are conjugates and others are not means
% that the imaginary part will always be the same but the half-life will
% differ between the two, which does not work out. To average intelligently
% pairs are detected, and then the imaginary and real part are averaged
% jointly, then the conjugates are re-formed as pairs. 

%% Average eigenvalues.

% Assign the average array.
eigValsAv = zeros(size(eigVals,1),1);

% Assign the currently examined eigenvalue.
examineEigVal = 1;

% Average over values.
while examineEigVal < size(eigVals,1)
    condition = sum(imag(eigVals(examineEigVal,:)) > 0) > 0;
    if condition
        if strcmp(flag,'mean')
            realPart = mean2(real(eigVals(examineEigVal:examineEigVal+1,:)));
            imaginaryPart = mean2(abs(imag(eigVals(examineEigVal:examineEigVal+1,:))));
        elseif strcmp(flag,'median')
            realPart = ...
                median(median(real(eigVals(examineEigVal:examineEigVal+1,:))));
            imaginaryPart = ...
                median(median(abs(imag(eigVals(examineEigVal:examineEigVal+1,:)))));
        end
        eigValsAv(examineEigVal) = complex(realPart,imaginaryPart);
        eigValsAv(examineEigVal+1) = complex(realPart,-imaginaryPart);
        examineEigVal = examineEigVal + 2;
    else
        if strcmp(flag,'mean')
            eigValsAv(examineEigVal) = mean2(real(eigVals(examineEigVal,:)));
        elseif strcmp(flag,'median')
            eigValsAv(examineEigVal) = median(real(eigVals(examineEigVal,:)));
        end
        examineEigVal = examineEigVal + 1;
    end
end

end