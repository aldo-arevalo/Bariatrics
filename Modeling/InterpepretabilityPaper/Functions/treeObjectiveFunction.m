function [objective] = treeObjectiveFunction(param, predictors,response,Options,cvp)

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationTree = fitctree(...
    predictors, ...
    response, ...
    'SplitCriterion', 'gdi', ...
    'Surrogate', 'on', ...
    'CrossVal','on',...
    'CVPartition',cvp, ...
    'Cost', [0 Options.Cost_FP; Options.Cost_FN 0], ...
    'ClassNames', [0; 1], ...
    'Prior', 'empirical', ...
    'MinLeafSize', param.MinLeafSize, ... 
    'MaxNumSplits', param.MaxNumSplits, ...
    'PredictorNames', Options.predictorNames, ...
    'ResponseName', Options.metric, ...
    'AlgorithmForCategorical', 'Exact', ...
    'CategoricalPredictors', Options.isCategoricalPredictor);

% Perform cross-validation
%partitionedModel = crossval(classificationTree, 'KFold', Options.k);

% Compute validation predictions
%[validationPredictions, validationScores] = kfoldPredict(partitionedModel);
[validationPredictions, validationScores] = kfoldPredict(classificationTree);

% Compute AUROC
[~,~,Ttree,~] = perfcurve(response,validationScores(:,end),1);
% [threshold , AUROC ] = getPerformanceMetrics(response,validationPredictions,Ttree);
[~ , AUROC ] = getPerformanceMetrics(response,validationPredictions,Ttree);

if Options.verbose
    [ACC,pr,speci,sens,Fs,~,] = confusion_matrix_([response validationPredictions]);
%     fprintf('Validation results Before adjusting threshold \n')
    fprintf('Accuracy %0.2f \n',ACC);
    fprintf('AUROC %0.2f \n', AUROC);
    fprintf('Precision %0.2f \n', pr);
    fprintf('Specificity %0.2f \n', speci);
    fprintf('Sensitivity %0.2f \n', sens);
    fprintf('F1 score %0.2f \n', Fs);
else
    fprintf('Accuracy %0.2f \n',ACC);
    fprintf('Sensitivity %0.2f \n', sens);
    fprintf('AUROC %0.2f \n', AUROC);
end

% compare_class = [response validationScores(:,end)];
% compare_class(compare_class(:, 2) >= threshold, 2) = 1;
% compare_class(compare_class(:, 2) < threshold, 2) = 0;
% [ACC,pr,speci,sens,Fs,~,] = confusion_matrix_(compare_class);
% [~,~,~,AUROC] = perfcurve(response,compare_class(:, 2),1);
% 
% fprintf('Validation results After adjusting threshold \n')
% fprintf('Accuracy %0.2f \n',ACC);
% fprintf('AUROC %0.2f \n', AUROC);
% fprintf('Precision %0.2f \n', pr);
% fprintf('Specificity %0.2f \n', speci);
% fprintf('Sensitivity %0.2f \n', sens);
% fprintf('F1 score %0.2f \n', Fs);

objective = AUROC;