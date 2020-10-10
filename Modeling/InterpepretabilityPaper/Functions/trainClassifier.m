function [trainedClassifier, partitionedModel, ...
    performance_results] = trainClassifier(trainingData, Options, ...
    MinLeafSize, MaxNumSplits, cvp)
% Returns trained classifiers and their performance.
%
%  Input:
%      trainingData: A table containing the same predictor and response
%       columns from the training/validation subset.
%
%  Output:
%      trainedClassifier: A struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier. This tree model can be used for making new predictions
%       on new data
%
%      trainedClassifier.predictFcn: A function to make predictions on new
%       data.
%
%      partitionModel: An object containing a
%       ClassificationPartitionedModel. This class contains a kFold cross 
%       validated model. Estimate the quality of classification by cross 
%       validation using one or more “kfold” methods: kfoldPredict, 
%       kfoldLoss, kfoldMargin, kfoldEdge, and kfoldfun. 
%       
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a classifier trained with the original data set
% T, enter:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 must be a table containing at least the same predictor columns as used
% during training. For details, enter:
%   trainedClassifier.HowToPredict

% Autor: Aldo Arévalo


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictors = inputTable(:, Options.predictorNames);
response = inputTable.(Options.metric);

%% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationTree = fitctree(...
    predictors, ...
    response, ...
    'SplitCriterion', 'gdi', ...
    'AlgorithmForCategorical', 'Exact', ...
    'CategoricalPredictors', Options.isCategoricalPredictor, ...
    'PredictorNames', Options.predictorNames, ...
    'ResponseName',Options.metric, ...
    'Surrogate', 'on', ...
    'Cost', [0 Options.Cost_FP; Options.Cost_FN 0], ...
    'Prior', 'empirical', ...
    'MinLeafSize', MinLeafSize, ... 
    'MaxNumSplits', MaxNumSplits, ...
    'ClassNames', [0; 1]);

%% Create the result struct with predict function

predictorExtractionFcn = @(t) t(:, Options.predictorNames);
treePredictFcn = @(x) predict(classificationTree, x);
trainedClassifier.predictFcn = @(x) treePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationTree = classificationTree;
trainedClassifier.About = 'This struct is a trained model.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

%% Perform cross-validation
% Training/validation for training subset

% This code specifies all the classifier options and trains the classifier
% kFold cross validation. This model can't be used for the predictions in
% the testData
partitionedModel = fitctree(...
    predictors, ...
    response, ...
    'SplitCriterion', 'gdi', ...
    'AlgorithmForCategorical', 'Exact', ...
    'CategoricalPredictors', Options.isCategoricalPredictor, ...
    'PredictorNames', Options.predictorNames, ...
    'ResponseName',Options.metric, ...
    'Surrogate', 'on', ...
    'CrossVal','on', ...
    'CVPartition', cvp,...
    'Cost', [0 Options.Cost_FP; Options.Cost_FN 0], ...
    'Prior', 'empirical', ...
    'MinLeafSize', MinLeafSize, ... 
    'MaxNumSplits', MaxNumSplits, ...
    'ClassNames', [0; 1]);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Performance metrics
% Compute performance metrics during training/validation
[AUROCval] = kfoldfun(partitionedModel, @treePerformanceKFold);

% Performance results struct
performance_results = [];
performance_results.AUROCval = AUROCval;
performance_results.meanAUROCval = mean(AUROCval,'omitnan');
performance_results.BestFold = find(AUROCval == max(AUROCval));
performance_results.ValidationPredictions = validationPredictions;
performance_results.ValidationPositive = validationScores;
performance_results.Mdl = classificationTree;

% Compute AUROC
[~,~,Ttree,~] = perfcurve(response,validationScores(:,end),1);
[threshold , ~ ] = getPerformanceMetrics(response,validationPredictions,Ttree);

performance_results.decisionBoundEst = threshold;

%% Save Tree view of each fold 
%Create folder Paritions in case it has not been created
foldername = 'Plots';
foldername2 = 'DecisionTree';
tmp_folder_exist = exist(foldername,'file');

if ne(tmp_folder_exist, 7)
    mkdir 'Plots'
end

tmp_folder_exist = exist([foldername,'/',foldername2],'file');

if ne(tmp_folder_exist, 7)
    mkdir Plots DecisionTree 
end

clearvars tmp_folder_exist

close all
delete(allchild(groot))

% Save tree view in classificationTree
view(classificationTree, 'Mode','Graph');
h = findall(groot,'Type','figure'); % Find all figures
filename = ['Tree_holdout_', Options.Patients,'_', num2str(Options.year),'yr'];
saveas(h,[foldername,'/',foldername2,'/',filename],'mfig');
saveas(h,[foldername,'/',foldername2,'/',filename],'eps');
delete(allchild(groot))

clearvars h

% Save each tree view in partitionedModel
for i=1:Options.k
    view(partitionedModel.Trained{i, 1}, 'Mode', 'Graph');
    h = findall(groot,'Type','figure'); % Find all figures
    filename = sprintf('Tree_%s_fold_%d_%d_yr',Options.Patients,i,Options.year);
    saveas(h,[foldername,'/',foldername2,'/',filename],'mfig');
    saveas(h,[foldername,'/',foldername2,'/',filename],'eps');
    delete(allchild(groot))
end

