function [resultsTS,selected_sets_TS] = featureSelectionFS(Options,trainingData,testData,cvPartition)

% Function to run Sequential Forward Selection (SFS) for Logistic
% Regression models
%
% INPUTS:
%
% OUTPUTS:
%
%
% Dependencies: M1_Modeling
%
% Author: Aldo Arévalo
% Date: 10/11/2020

% Get names for output variable and input variables
responseName = Options.metric;
predictorNames = Options.predictorNames;

% Remove categorical dtype
trainingData = convertvars(trainingData,@iscategorical,'double');
testData = convertvars(testData,@iscategorical,'double');

% Test subset
Xtest = testData(:,predictorNames);
Ytest = testData(:,responseName);

% CV
kfolds = Options.k;

choice = Options.SFSOptions;

count=1;

for ii=1:kfolds
    
    cv_train = trainingData(training(cvPartition, ii),:); % input data for training
    cv_val = trainingData(test(cvPartition, ii),:); % input data for validation
 
    %% SFS with TakagiSugeno-FCM
    
    if Options.Fuzzy
        % Define also the structure with model parameters
        FM.c = Options.cluster;        % number of clusters
        FM.m =  2;       % fuzziness parameter
        FM.tol = 0.01;   % termination criterion
        FM.ante = 2;     % type of antecedent:   1 - product-space MFS 2 - projected MFS
        FM.cons = 2;
        FM.seed = Options.seed;
        
        % Training subset
        Xtrain = cv_train(:,predictorNames);
        Ytrain = cv_train(:,responseName);
        
        Xval = cv_val(:,predictorNames);
        Yval = cv_val(:,responseName);
        
%         try
            [selected_set, AUC_selected_set, results] = ...
                SFSTS(Xtrain,Ytrain,Xval,Yval,choice,FM,predictorNames,...
                responseName,Options.clust_method,Options.seed, ...
                Options.undersampling_factor);
            
            if choice(3)==1
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                disp(['Takagi-Sugeno: '])
                disp(predictorNames(selected_set))
                disp(['AUC: ', num2str(round(AUC_selected_set,2))])
                disp(['Number of iterations/models created: ', num2str(size(results,1))])
            end
            
            close all
            
            selected_sets_TS(count)={selected_set};
            
            % TRAIN model
            Dat.U=table2array(trainingData(:,selected_set));
            Dat.Y=table2array(trainingData(:,responseName));
            Dat.InputName = predictorNames(selected_set);
            Dat.OutputName = responseName;
            
            if strcmp(Options.clust_method,'FGK')
                [model, ~] = fmclust_silent(Dat,FM);
            end

            if strcmp(Options.clust_method,'FCM')
                [model, ~] = fmclust2(Dat,FM);
            end
            
            yfitTrain = fmsim(Dat.U,Dat.Y,model,[],[],0);
            
            % Optimize the threshold
            [~,~,Thresholds,~] = perfcurve(table2array(...
                trainingData(:,responseName)),...
                yfitTrain,1);
            [threshold , ~ ] = getPerformanceMetrics(...
                table2array(trainingData(:,responseName)),...
                yfitTrain,Thresholds);
            
            % TESTING model
            Dat.U=table2array(Xtest(:,selected_set));
            Dat.Y=table2array(Ytest);
            
            yfitTest = fmsim(Dat.U,Dat.Y,model,[],[],0);
            
            [X,Y,~,AUROC] = perfcurve(testData.TWL, yfitTest, 1);
        
            % Plot AUROC
            plot(X,Y)
            
            compare_class = [Ytest.TWL yfitTest];
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
            
            resultsTS(count).set = selected_set;
            resultsTS(count).VarNames = predictorNames(selected_set);
            resultsTS(count).AUC = AUROC;
            resultsTS(count).accuracy = ACC;
            resultsTS(count).sensitivity = sens;
            resultsTS(count).specificity = speci;
            resultsTS(count).F1 = Fs;
            resultsTS(count).precision = pr;
            resultsTS(count).threshold=threshold;
            resultsTS(count).ConfusionMat = cfxmat;
            resultsTS(count).Model = model;
            
%         catch
%             fismat = [];
%             
%             resultsTS(count).set = [];
%             resultsTS(count).VarNames = {};
%             resultsTS(count).AUC = NaN;
%             resultsTS(count).accuracy = NaN;
%             resultsTS(count).sensitivity = NaN;
%             resultsTS(count).specificity = NaN;
%             resultsTS(count).threshold = NaN;
%             resultsTS(count).F1 = NaN;
%             resultsTS(count).precision = NaN;
%             resultsTS(count).FM = FM;
%             resultsTS(count).ConfusionMat = NaN;
%         end
    end

    count=count+1;
    close all
end

% The function returns:
% selected_set
% AUC_selected_set
% results: a matrix with all the feature combinations tested during SFS, from the
% first iteration to the last. The last column gives the AUC for each
% combination.

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
