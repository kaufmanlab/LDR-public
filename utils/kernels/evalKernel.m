function kernelMat = evalKernel(measurements1,measurements2,type,lengthscale)
%% OVERVIEW

% This function parameterizes the similarity measure between two
% measurements.

%% EVALUATE KERNEL.

% Evaluate the kernel for distance kernels.
if strcmp(type,'squaredExponential')
    kernelMat = pdist2(measurements1.',measurements2.')/lengthscale;
    kernelMat = exp(-kernelMat.^2);
elseif strcmp(type,'OU')
    kernelMat = pdist2(measurements1.',measurements2.')/lengthscale;
    kernelMat = exp(-abs(kernelMat));
end

end

