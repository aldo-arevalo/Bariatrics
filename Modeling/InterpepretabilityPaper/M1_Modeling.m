%% Main modeling

clearvars
close all

addpath(genpath('Datasets'))
addpath(genpath('Functions'))

% Options
Kappa = 0;
AUC = 1;
Reg = 0;
bala = 1; % True only for Classification models

% Control options
prova = Kappa + AUC + Reg;
if prova>=2
    error('Error. Select only one option, not 2 or more.')
end

%% Fuzzy Modeling
% filePatternFuzzy = fullfile(foldername1,foldername2,foldername3,{'FCM*.mat';'FGK*.mat'});
% 
% for p=1:size(filePatternFuzzy,1)
%     matFilesFuzzy = dir(filePatternFuzzy{p,1});
%     for ko=1:length(matFilesFuzzy)
%         matFilename = fullfile(foldername1,foldername2,foldername3, matFilesFuzzy(ko).name);
%         load(matFilename,'Options','z','text','metric','X','Y','results',...
%             'labels','filename_partition');
%     close all
%         M5_ModelingReg
%         clearvars response
%     end
% end

%% Logistic Regression

% filePatternLR = fullfile(foldername1,foldername2,foldername3,{'LR*.mat'});
% 
% for p=1:size(filePatternLR,1)
%     matFilesLR = dir(filePatternLR{p,1});
%     for ko=1:length(matFilesLR)
%         matFilename = fullfile(foldername1,foldername2,foldername3, matFilesLR(ko).name);
%         load(matFilename,'Options','z','text','metric','X','Y','results',...
%             'labels','filename_partition');
%     close all
%         M5_ModelingReg
%         clearvars response
%     end
% end

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

%% Decision Tree
tic
% No SFS is needed for DT

surgery_type = {'All','Metab','Bar'};

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = false;
Options.Fuzzy = false;
Options.SVM = false;
Options.DT = true;
Options.Kappa = false;
Options.algth = 'DT';
% CV folds
Options.k = 5;
% Hold out
Options.Hold_out = 0.2;
% Partitions (runs)
Options.n_partitions = 'No SFS needed';
% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'
% Other Options
Options.cluster = 'No apply';
Options.clust_method = 'No apply';
Options.Regression = false;
Options.Classification = true;
Options.JOB = 'Interpret';
Options.seed = 13;
Options.Cost_FP = 5;
Options.Cost_FN = 1;
Options.verbose = true;

filename = '(2020-10-08 06 41 29) CZE(1yr&2yr)_25%TWL_53var_.mat';

rng(Options.seed); % For reproducibility

to_save = 'FilesModeled_DT';

% Building strings
foldername1 = 'Results';

% Create folder Results in case it has not been created
tmp_folder_exist = exist(foldername1,'file');
    
if ne(tmp_folder_exist, 7)
    mkdir 'Results'
end

clearvars tmp_folder_exist

foldername2 = 'AUC';
% Create folder AUC in case it has not been created
tmp_folder_exist = exist([foldername1,'/',foldername2],'file');

if ne(tmp_folder_exist, 7)
    mkdir 'Results/AUC'
end

clearvars tmp_folder_exist

for i=1:length(surgery_type)
    Options.Patients = surgery_type{i};

    for a=1:2
        year = a;

        switch year
            case 1
                if strcmp(Options.Patients,'All')
                    load(filename,'tbl_in_1yr','var_names1')
                    Options.input = tbl_in_1yr;
                end
                if strcmp(Options.Patients,'Metab')
                    load(filename,'tbl_1yr_metabolic','var_names1')
                    Options.input = tbl_1yr_metabolic;
                end
                if strcmp(Options.Patients, 'Bar')
                    load(filename,'tbl_1yr_bariatric','var_names1')
                    Options.input = tbl_1yr_bariatric;
                end
                
                Options.year = year;
                idx_patient = strcmp(var_names1, 'PatientCode');
                var_names1(idx_patient) = [];
                idx_patient = strcmp(var_names1, 'dyslip');
                var_names1(idx_patient) = [];
                idx_patient = strcmp(var_names1, 'hypert');
                var_names1(idx_patient) = [];
                idx_patient = strcmp(var_names1, 'Age_');
                var_names1(idx_patient) = [];
                Options.var_names_original = var_names1;
                
                clearvars idx_patient var_names1
                
            case 2
                if strcmp(Options.Patients,'All')
                    fprintf('All patients')
                    load(filename,'tbl_in_2yr','var_names2')
                    Options.input = tbl_in_2yr;
                end
                if strcmp(Options.Patients,'Metab')
                    fprintf('Metabolic patients')
                    load(filename,'tbl_2yr_metabolic','var_names2')
                    Options.input = tbl_2yr_metabolic;
                end
                if strcmp(Options.Patients, 'Bar')
                    fprintf('Bariatric patients')
                    load(filename,'tbl_2yr_bariatric','var_names2')
                    Options.input = tbl_2yr_bariatric;
                end
                Options.year = year;
                idx_patient = strcmp(var_names2, 'PatientCode');
                var_names2(idx_patient) = [];
                idx_patient = strcmp(var_names2, 'dyslip');
                var_names2(idx_patient) = [];
                idx_patient = strcmp(var_names2, 'hypert');
                var_names2(idx_patient) = [];
                idx_patient = strcmp(var_names2, 'Age_');
                var_names2(idx_patient) = [];
                Options.var_names = var_names2;
                
                clearvars var_names2 idx_patients
        end
        
        % Prepare data for modeling
        [holdPartition,trainingData,testData,cvPartition,Options] = getDatasetTree(Options);
        
        %% Train and optimize parameters
        
        % Maximize AUROC
        [BayesResults] = treeBayesOptimize(trainingData,cvPartition,Options);
        close all
        
        % Extract best parameters
        MinLeafSize = BayesResults.XAtMinEstimatedObjective.MinLeafSize;
        MaxNumSplits = BayesResults.XAtMinEstimatedObjective.MaxNumSplits;
        
        % Train model
        [trainedClassifier, partitionedModel, performance_results] = trainClassifier(...
            trainingData, Options, MinLeafSize, MaxNumSplits, cvPartition);
        
        %% Test model

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
        
        file_to_save = [filename(1:48),Options.Patients,'_', ...
            num2str(Options.year),'yr','.mat'];
        save(fullfile(foldername1,foldername2,file_to_save), ...
            'performance_results', 'trainedClassifier', ...
            'partitionedModel', 'BayesResults', 'MaxNumSplits', ...
            'MaxNumSplits', 'Options')
        
        clearvars pr sens speci ACC AUC AUROC Fs

        close all
    end
end

fprintf('Decision Tree models finished \n')

clearvars i

% Save the files already modeled

foldername3 = 'Workspaces';

% Create folder Workspaces in case it has not been created
tmp_folder_exist = exist([foldername1,'/',foldername2,'/', ...
    foldername3],'file');

if ne(tmp_folder_exist, 7)
    mkdir 'Results/AUC/Workspaces'
end

clearvars tmp_folder_exist

fprintf('Elapsed time is %0.2f minutes.',toc/60)