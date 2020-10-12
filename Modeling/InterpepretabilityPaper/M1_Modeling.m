%% Main modeling

clearvars
close all

addpath(genpath('Datasets'))
addpath(genpath('Functions'))

%% Define general setting for the models
% Requires manual input

% File that contains the datasets
filename = '(2020-10-08 06 41 29) CZE(1yr&2yr)_25%TWL_53var_.mat';

% For reproducibility
Options.seed = 13;
% CV folds
Options.k = 5;
% Hold out
Options.Hold_out = 0.2;
% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

% Clusters for Fuzzy models
clusters = [2 , 3];
clustering = {'FGK','FCM'};
% FGK - Gustaffson-Kessel
% FCM - Fuzzy C-means

% Objective function
Kappa = 0;
AUC = 1;
Fscore = 0;

% Define SFS options
choice(1)=1; % stop criterion:
% if 0 % dont stop when there is no improvement
% if 1 % stop when there is no improvement (faster)

choice(2)=1; % plot results:
% if 0 % dont plot
% if 1 % plot performance results and number of features

choice(3)=1; % display the results
% if 0 % dont display
% if 1 % display

choice(4)=1; % random undersampling without replacement
% if 0 % no class undersampling
% if 1 % class undersampling
% Select a factor for undersampling
% EXAMPLES:
%   1 => 1 entry from minority class:1 entry from dominant class
%   1.3 => 1 entry from minority class:1.3 entries from dominant class
Options.undersampling_factor = 1.3;

% Penalization for false positives (FP) and false negatives (FN)
Options.Cost_FP = 5;
Options.Cost_FN = 1;

% Other options
Options.Regression = false;
Options.Classification = true;
Options.JOB = 'Interpret';

% Display results for Bayesian Optimization
Options.verbose = true;

%% Automatic settings
prova = Kappa + AUC + Fscore;
if prova>=2
    error('Error. Select only one Objective Function, not 2 or more.')
end

% For reproducibility
rng(Options.seed,'twister'); 

% Surgery groups
% 3 groups of patients modeled
surgery_type = {'All','Metab','Bar'};

% Building strings
foldername1 = 'Results';

% Create folder Results in case it has not been created
tmp_folder_exist = exist(foldername1,'file');
    
if ne(tmp_folder_exist, 7)
    mkdir 'Results'
end

clearvars tmp_folder_exist prova

% Objective function, what is optimized?
if AUC == 1
    foldername2 = 'AUC';
    % Create folder AUC in case it has not been created
    tmp_folder_exist = exist([foldername1,'/',foldername2],'file');
    if ne(tmp_folder_exist, 7)
        if choice(4) == 1
            foldername3 = 'Undersampling';
            % Create folder Undersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/AUC/Undersampling'
            end
        else
            foldername3 = 'NoUndersampling';
            % Create folder NoUndersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/AUC/NoUndersampling'
            end
        end
    else
        if choice(4) == 1
            foldername3 = 'Undersampling';
            % Create folder Undersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/AUC/Undersampling'
            end
        else
            foldername3 = 'NoUndersampling';
            % Create folder NoUndersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/AUC/NoUndersampling'
            end
        end
    end
end

if Fscore == 1
    foldername2 = 'Fscore';
    % Create folder Fscore in case it has not been created
    tmp_folder_exist = exist([foldername1,'/',foldername2],'file');
    if ne(tmp_folder_exist, 7)
        if choice(4) == 1
            foldername3 = 'Undersampling';
            % Create folder Undersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/Fscore/Undersampling'
            end
        else
            foldername3 = 'NoUndersampling';
            % Create folder NoUndersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/Fscore/NoUndersampling'
            end
        end
    else
        if choice(4) == 1
            foldername3 = 'Undersampling';
            % Create folder Undersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/Fscore/Undersampling'
            end
        else
            foldername3 = 'NoUndersampling';
            % Create folder NoUndersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/Fscore/NoUndersampling'
            end
        end
    end
end

if Kappa == 1
    foldername2 = 'Kappa';
    % Create folder Kappa in case it has not been created
    tmp_folder_exist = exist([foldername1,'/',foldername2],'file');
    if ne(tmp_folder_exist, 7)
        if choice(4) == 1
            foldername3 = 'Undersampling';
            % Create folder Undersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/Kappa/Undersampling'
            end
        else
            foldername3 = 'NoUndersampling';
            % Create folder NoUndersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/Kappa/NoUndersampling'
            end
        end
    else
        if choice(4) == 1
            foldername3 = 'Undersampling';
            % Create folder Undersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/Kappa/Undersampling'
            end
        else
            foldername3 = 'NoUndersampling';
            % Create folder NoUndersampling in case it has not been created
            tmp_folder_exist2 = exist([foldername1,'/',foldername2,'/',foldername3],'file');
            if ne(tmp_folder_exist2, 7)
                mkdir 'Results/Kappa/NoUndersampling'
            end
        end
    end
end

Options.ObjFunction = foldername2;

clearvars tmp_folder_exist AUC Fscore Kappa tmp_folder_exist2

%% SVM

% filePatternSVM = fullfile(foldername1,foldername2,foldername3,{'SVM*.mat'});
% 
% for p=1:size(filePatternSVM,1)
%     matFilesSVM = dir(filePatternSVM{p,1});
%     for ko=1:length(matFilesSVM)
%         matFilename = fullfile(foldername1,foldername2,foldername3, matFilesSVM(ko).name);
%         load(matFilename,'Options','z','text','metric','X','Y','results',...
%             'labels','filename_partition');
%     close all
%         M5_ModelingSVM
%         clearvars response
%     end
% end

%% Models
tic

for i=1:length(surgery_type)
    Options.Patients = surgery_type{i};

    for a=1:2
        Options.year = a;
        
        %% Decision Trees
        % SFS is not needed for DT
        Options.SFSOptions = 'SFS is not needed for DT';
        
        %> Select Model options
        Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
        Options.Logistic = false;
        Options.Fuzzy = false;
        Options.SVM = false;
        Options.DT = true;
        Options.algth = 'DT';

        %> Other Options
        Options.cluster = 'No apply';
        Options.clust_method = 'No apply';
        
        [holdPartition,trainingData,testData,...
            cvPartition,Options] =  getInputSubsets(Options, filename);
        
        %> Train and optimize parameters
            % Maximize AUROC
            [BayesResults] = treeBayesOptimize(trainingData,cvPartition,Options);
            close all

            % Extract best parameters
            MinLeafSize = BayesResults.XAtMinEstimatedObjective.MinLeafSize;
            MaxNumSplits = BayesResults.XAtMinEstimatedObjective.MaxNumSplits;
        
        %> Train model
            [trainedClassifier, partitionedModel, performance_results] = trainClassifier(...
                trainingData, Options, MinLeafSize, MaxNumSplits, cvPartition);
        
        %> Test model

            % Compute predictions for test set
            [yfit_test, yscores_test] = trainedClassifier.predictFcn(...
                testData);

            [X,Y,Ttree,AUROC] = perfcurve(testData.TWL, yscores_test(:, 2),1);
            performance_results.AUROCbefore = AUROC;
        
            % Plot AUROC
            plot(X,Y)

            [ACC,pr,speci,sens,Fs,~,] = confusion_matrix_([testData.TWL yfit_test]);
            fprintf('Validation results Before adjusting threshold \n')
            fprintf('Accuracy %0.2f \n',ACC); performance_results.ACCbefore = ACC;
            fprintf('AUROC %0.2f \n', AUROC); 
            fprintf('Precision %0.2f \n', pr); performance_results.PrecisionBefore = pr;
            fprintf('Specificity %0.2f \n', speci); performance_results.SpecifBefore = speci;
            fprintf('Sensitivity %0.2f \n', sens); performance_results.SensiBefore = sens;
            fprintf('F1 score %0.2f \n', Fs); performance_results.F1Before = Fs;

            % Confusion matrix
            figure
            cm = confusionchart(testData.TWL, yfit_test);
            cm.Title = 'Weight Loss Prediction';
            cm.RowSummary = 'row-normalized';
            cm.ColumnSummary = 'column-normalized';

            threshold = performance_results.decisionBoundEst;
            compare_class = [testData.TWL yscores_test(:, 2)];
            compare_class(compare_class(:, 2) >= threshold, 2) = 1;
            compare_class(compare_class(:, 2) < threshold, 2) = 0;
            [ACC,pr,speci,sens,Fs,~,] = confusion_matrix_(compare_class);
            [~,~,~,AUC] = perfcurve(testData.TWL,compare_class(:, 2),1);
            performance_results.ValidationPositiveThreshold = compare_class(:, 2);

            fprintf('Validation results After adjusting threshold \n')
            fprintf('Accuracy %0.2f \n',ACC); performance_results.ACCafter = ACC;
            fprintf('AUROC %0.2f \n', AUC); performance_results.AUROCafter = AUC; 
            fprintf('Precision %0.2f \n', pr); performance_results.PrecisionAfter = pr;
            fprintf('Specificity %0.2f \n', speci); performance_results.SpecifAfter = speci;
            fprintf('Sensitivity %0.2f \n', sens); performance_results.SensiAfter = sens;
            fprintf('F1 score %0.2f \n', Fs); performance_results.F1After = Fs;

            % Confusion matrix
            figure
            cm = confusionchart(testData.TWL,compare_class(:, 2));
            cm.Title = 'Weight Loss Prediction Afteradjusting Threshold';
            cm.RowSummary = 'row-normalized';
            cm.ColumnSummary = 'column-normalized';

            file_to_save = [filename(1:22),Options.algth,'_',...
                Options.Patients,'_', num2str(Options.year),'yr','.mat'];
            save(fullfile(foldername1,foldername2,foldername3,file_to_save), ...
                'performance_results', 'trainedClassifier', ...
                'partitionedModel', 'BayesResults', 'MaxNumSplits', ...
                'MaxNumSplits', 'Options','cvPartition','holdPartition')

            clearvars pr sens speci ACC AUC AUROC Fs

        close all
        fprintf('Decision Tree models finished \n')
        fprintf('Elapsed time is %0.2f minutes.',toc/60)
        clearvars BayesResults cm compare_class file_to_save MaxNumSplits ...
            MinLeafSize performance_results threshold trainedClassifier ...
            Ttree X Y yfit_test yscores_test partitionedModel
        %% Logistic regression
        
        %> Select Model options
        Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
        Options.Logistic = true;
        Options.Fuzzy = false;
        Options.SVM = false;
        Options.DT = false;
        Options.algth = 'LR';
        
        % SFS options
        Options.SFSOptions = choice;
        
        % Run SFS
        [resultsLR, ~] = featureSelectionLR(Options, trainingData, testData, cvPartition);
        
        file_to_save = [filename(1:22),Options.algth,'_',...
                Options.Patients,'_', num2str(Options.year),'yr','.mat'];
        save(fullfile(foldername1,foldername2,foldername3,file_to_save), ...
            'Options','resultsLR','cvPartition','holdPartition')
        
        clearvars resultsLR
        %% Fuzzy modeling
        
        %> Select Model options
        Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
        Options.Logistic = true;
        Options.Fuzzy = true;
        Options.SVM = false;
        Options.DT = false;
        Options.algth = 'F_TS';
        
        for c=clusters
            % Number of clusters
            Options.cluster = c;
            
            % Clustering algorithm
            for cc=1:length(clustering)
                Options.clust_method = clustering{cc};
                
                [resultsTS, ~] = featureSelectionFS(Options, trainingData, testData, cvPartition);
                
                file_to_save = [filename(1:22),...
                    Options.clust_method,'_',num2str(c),'c_',...
                    Options.Patients,'_', num2str(Options.year),...
                    'yr','.mat'];
                save(fullfile(foldername1,foldername2,foldername3,file_to_save), ...
                    'Options','resultsTS','cvPartition','holdPartition')
            end
        end
        
        fprintf('Fuzzy models finished \n')
        fprintf('Elapsed time is %0.2f minutes.',toc/60)
        
    end
end

clearvars i

fprintf('Modeling completed \n')
fprintf('Elapsed time is %0.2f minutes.',toc/60)
%% Save 
% Save the files already modeled

foldername4 = 'Workspaces';

% % Create folder Workspaces in case it has not been created
% tmp_folder_exist = exist([foldername1,'/',foldername2,'/', foldername3,'/',...
%     foldername4],'file');
% 
% if ne(tmp_folder_exist, 7)
%     mkdir 'Results/AUC/Workspaces'
% end
% 
% clearvars tmp_folder_exist

fprintf('Elapsed time is %0.2f minutes.',toc/60)