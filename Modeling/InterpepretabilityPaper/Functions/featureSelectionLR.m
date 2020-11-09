function [resultsLR,selected_sets_LR] = featureSelectionLR(Options,trainingData,testData,cvPartition)

% Function to run Sequential Forward Selection (SFS) for Logistic
% Regression models
%
% INPUTS:
%       Options:      Structure that contains the general settings.
%       trainingData: Table containing the input subset for training.
%       testData:     Table containing the input subset for testing.
%       cvPartition:  Stratified cross validation partition.
%
% OUTPUTS:
%       resultsLR: Structure (1 x kFolds) with 10 fields summarizing the
%       indices of the selected features in each fold, their respective
%       variable names, the performance metrics in each fold, the decision
%       threshold to categorize the estimated output, the confusion
%       matrix, and the GLM object.
%       selected_sets_LR: Cell array  with the size {1 x kFold} containing
%                         the indices for the selected feature subseat in 
%                         each fold.
%
%
% Dependencies: M1_Modeling
%
% Author: Aldo Arévalo
% Date: 10/11/2020

% Get names for output variable and input variables
responseName = Options.metric;
predictorNames = Options.predictorNames;

% CV
kfolds = Options.k;

choice = Options.SFSOptions;

count=1;

for ii=1:kfolds
    
    cv_train = trainingData(training(cvPartition, ii),:); % input data for training
    cv_val = trainingData(test(cvPartition, ii),:); % input data for validation
    
    %% SFS with Logistic Regression
    if Options.Logistic
        %> SFS with LR

        [selected_set, AUC_selected_set, results] = ...
            SFSLR(cv_train,cv_val,choice,responseName,predictorNames,...
            Options.undersampling_factor);

        if Options.verbose
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(['Logistic regression: '])
            disp(predictorNames(selected_set))
            disp(['AUC: ', num2str(round(AUC_selected_set,2))])
            disp(['Number of iterations/models created: ', num2str(size(results,1))])
        end
        
        close all
        
        selected_sets_LR(count) = {selected_set};
        var_names = predictorNames(selected_set);
        var_names{end+1} = responseName;

        % TRAIN complete training dataset
        model = fitglm(trainingData(:,var_names), ...
            'Distribution', 'binomial', ...
            'link', 'logit',...
            'ResponseVar',responseName);
                    
        yfitTrain = predict(model,trainingData(:,var_names));
        
        % Optimize threshold
        [~,~,Thresholds,~] = perfcurve(trainingData.(responseName),...
            yfitTrain,1);
        [threshold , ~ ] = getPerformanceMetrics(trainingData.(responseName),...
            yfitTrain,Thresholds);

        % FINAL test
        yscoresTest = predict(model,testData(:,var_names));
        
        [X,Y,~,AUROC] = perfcurve(testData.TWL, yscoresTest, 1);
        
        % Plot AUROC
        plot(X,Y)
        
        compare_class = [testData.TWL yscoresTest];
        compare_class(compare_class(:, 2) >= threshold, 2) = 1;
        compare_class(compare_class(:, 2) < threshold, 2) = 0;
        
        [ACC,pr,speci,sens,Fs,~,cfxmat] = confusion_matrix_(compare_class);
        fprintf('Validation results adjusting threshold \n')
        fprintf('Accuracy %0.2f \n',ACC);
        fprintf('AUROC %0.2f \n', AUROC); 
        fprintf('Precision %0.2f \n', pr);
        fprintf('Specificity %0.2f \n', speci);
        fprintf('Sensitivity %0.2f \n', sens);
        fprintf('F1 score %0.2f \n', Fs);
            
        resultsLR(count).set = selected_set;
        resultsLR(count).VarNames = predictorNames(selected_set);
        resultsLR(count).AUC = AUROC;
        resultsLR(count).accuracy = ACC;
        resultsLR(count).sensitivity = sens;
        resultsLR(count).specificity = speci;
        resultsLR(count).F1 = Fs;
        resultsLR(count).precision = pr;
        resultsLR(count).threshold=threshold;
        resultsLR(count).ConfusionMat = cfxmat;
        resultsLR(count).Model = model;
        
        close all
    end

    count=count+1;
end


%% plot SFS final results
% selected_LR=zeros(1,size(Xtrain,2));
% selected_logical_LR=zeros(kfolds,size(input,2));
% selected_TS=zeros(1,size(Xtrain,2));
% selected_logical_TS=zeros(kfolds*rep_cross,size(input,2));
% for i=1:kfolds%*rep_cross
%     aux_LR=cell2mat(selected_sets_LR(1,i));
%     %aux_TS=cell2mat(selected_sets_TS(1,i));
%     selected_logical_LR(i,aux_LR)=1;
%     %selected_logical_TS(i,aux_TS)=1;
%     for ii=1:size(Xtrain,2)
%         if isempty(find(aux_LR==ii))
%         else
%             selected_LR(ii)=selected_LR(ii)+1;
%         end
%         %if isempty(find(aux_TS==ii))
%         %else
%         %    selected_TS(ii)=selected_TS(ii)+1;
%         %end
%     end
% end
% 
% ftsize=14;
% plot_SFS(selected_logical_LR,predictorNames,ftsize);
% plot_SFS(selected_logical_TS,predictorNames,ftsize);
