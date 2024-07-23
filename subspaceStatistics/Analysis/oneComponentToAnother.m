function output = oneComponentToAnother
%% OVERVIEW

% This function attempts to predict one portion of the activity manifold of
% motor cortex from another. 
 
%% Parameters.

% Refers to the regularization on the regression.
L2 = 0.1;

% Refers to the portion of the singular value spectrum of the initial
% states to keep.
svEnergy = 0.8;

% Refers to the portion of the condition to leave out from
% cross-validation.
leaveOut = 6;

%% Fit dynamics.

% Load data.
load('ShenoyMonkeyData');
ShenoyMonkeyData = ShenoyMonkeyData(1:2);
load('basisFxnNum');

% Loop over monkeys. 
for monkey = 1:size(ShenoyMonkeyData,2)
    pooled = poolDatasets(ShenoyMonkeyData(monkey).M1,ShenoyMonkeyData(monkey).PMd);
    output(monkey).M1 = analyze(ShenoyMonkeyData(monkey).M1,L2,svEnergy, ...
        leaveOut,basisFxnNum(monkey).M1.maxDim);
    output(monkey).PMd = analyze(ShenoyMonkeyData(monkey).PMd,L2,svEnergy, ...
        leaveOut,basisFxnNum(monkey).PMd.maxDim);
end

end

%% SUBFUNCTION FOR ANALYZING DATA.

function out = analyze(data,L2,svEnergy,leaveOut,factorDim)
data = pruneRepeats(data);
% Perform prediction.
[out.data,out.prediction,reassambleTensor] = predictOneComponentFromAnother( ...
    data,leaveOut,L2,factorDim,svEnergy);
% Quantify the variance explained.
% for predictorIndex = 1:ceil(factorDim/2)
%     for predictedIndex = 1:ceil(factorDim/2)
%         for cond = 1:size(data,2)
%             out.varExplained(predictorIndex,predictedIndex,cond) = ...
%                 getAlignmentIndex(out.data(predictedIndex,cond).matrix, ...
%                 out.prediction(predictorIndex,predictedIndex,cond).matrix);
%         end
%     end
% end
for predictorIndex = 1:ceil(factorDim/2)
    for predictedIndex = 1:ceil(factorDim/2)
        out.varExplained(predictorIndex,predictedIndex) = ...
            getVarExplained(out.data(predictedIndex,:), ...
            squeeze(out.prediction(predictorIndex,predictedIndex,:)).','ind');
    end
end
% Shuffle the relationships between submanifolds.
out.nullData = data;
for cond = 1:size(out.nullData,2)
    out.nullData(cond).matrix = 0;
    for comp = 1:ceil(factorDim/2)
        randCond = randi([1 size(out.nullData,2)],1,1);
        out.nullData(cond).matrix = out.nullData(cond).matrix + ...
            reassambleTensor(comp,randCond).matrix;
    end
end
% Predict activity from shuffled initial states.
[out.nullData,out.nullPrediction] = predictOneComponentFromAnother( ...
    out.nullData,leaveOut,L2,factorDim,svEnergy);
% Quantify the variance explained.
% for predictorIndex = 1:ceil(factorDim/2)
%     for predictedIndex = 1:ceil(factorDim/2)
%         for cond = 1:size(data,2)
%             out.nullVarExplained(predictorIndex,predictedIndex,cond) = ...
%                 getAlignmentIndex(out.nullData(predictedIndex,cond).matrix, ...
%                 out.nullPrediction(predictorIndex,predictedIndex,cond).matrix);
%         end
%     end
% end
for predictorIndex = 1:ceil(factorDim/2)
    for predictedIndex = 1:ceil(factorDim/2)
        out.nullVarExplained(predictorIndex,predictedIndex) = ...
            getVarExplained(out.nullData(predictedIndex,:), ...
            squeeze(out.nullPrediction(predictorIndex,predictedIndex,:)).','ind');
    end
end
% Quantify the p value.
for predictorIndex = 1:ceil(factorDim/2)
    for predictedIndex = 1:ceil(factorDim/2)
        out.pVal(predictorIndex,predictedIndex) = ranksum(...
            out.nullVarExplained(predictorIndex,predictedIndex).array, ...
            out.varExplained(predictorIndex,predictedIndex).array);
    end
end
% Get summary.
out.meanArray = zeros(ceil(factorDim/2));
out.meanArrayNull = out.meanArray;
for predictorIndex = 1:ceil(factorDim/2)
    for predictedIndex = 1:ceil(factorDim/2)
        out.meanArray(predictorIndex,predictedIndex) = ...
            out.varExplained(predictorIndex,predictedIndex).mean;
        out.meanArrayNull(predictorIndex,predictedIndex) = ...
            out.nullVarExplained(predictorIndex,predictedIndex).mean;
    end
end
end

