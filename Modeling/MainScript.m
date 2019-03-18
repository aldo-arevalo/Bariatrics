%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dataset: Catharina Hospital
% Language: MATLAB 2016 a and beyond
% Description: Mainscript used to build Fuzzy models, Logistic Regression (LR) and SVM
% Input Variables:
%    Options                             Structure that contains all the input information needded to create the models
%                                        aforementioned. Customize the options according to the addressed model.
%    Options.JOB                         (str) Main folder where is going to be stored the files and plots
%    Options.Polynomial                  (logical) set as true if Univariable linear regression is dessired
%    Options.Logistic                    (Logical) set as true if LR is desired
%    Options.Fuzzy                       (Logical) set as true if Fuzzy modeling is desired.
%    Options.SVM                         (Logical) set as true if SVM model is desired
%    Options.Kappa                       (Logical) set as true if Cohen's coefficient is desired to estimate AUK as per-
%                                         formance metric
%    Options.Patients                    (str) 'All' all types of patiens; 'Metab' only metabolic patients; 'Bar' only 
%                                         bariatric patients
%    Options.Regression                  (logical) set as true if regression models are desired
%    Options.Classification              (logical) set as true if classification models are desired
%    Options.input                       (table) Should contain the input dataset including patient code column
%    Options.year                        (int) Should contain the follow-up year
%    Options.var_names                   (cell,str) Should contain the names of all columns of input.
%
% %  Options for k-fold cross validation
%    Options.n_partitions                (int) Set the number of runs or partitions needed. Each partition will be vali-
%                                         dated through k-fold cross validation
%    Options.k                           (int) Set the number of folds for cross validation
%
%    % Output metric                     Set only one as true
%    Options.TWL                         (logical) set as true if Total weight loss is desired as output metric
%    Options.EBWL                        (logical) set as true if Excess body weight loss is desired as output metric
%    Options.EBMIL                       (logical) set as true if Excess BMI loss is desired as output metric
%    Options.WL                          (logical) set as true if Absolute weight loss is desired as output metric
%
%    % Fuzzy modeling definitions
%    Options.clust_method                (str) 'FCM' Fuzzy c-means clustering
                                               'FGK' Gustaffson-Kessel clustering algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preamble
clearvars
addpath('SetnesII')
addpath('classification_toolbox')
addpath('fmid-v40')
addpath('Datasets/Deletion')

% Options is a structure that contains all the input information needed to create the models
Options.JOB = 'SetnesII';

%% Fuzzy Regression models - All patients

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = false;
Options.Fuzzy = true;
Options.SVM = false;
Options.Kappa = false;
Options.Patients = 'All';

% CV fold
Options.k = 5;

% Partitions (runs)
Options.n_partitions = 10;

% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

for a=1:2
    year = a;
    for c=2:6
        Options.cluster = c;
        Options.Regression = true;
        Options.Classification = false;
        
        switch year
            case 1
                load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
                %load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                    'tbl_in_1yr','var_names1')
                Options.input = tbl_in_1yr;
                Options.year = year;
                Options.var_names = var_names1;
            case 2
                load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
                %load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                      'tbl_in_2yr','var_names2')
                Options.input = tbl_in_2yr;
                Options.year = year;
                Options.var_names = var_names2;
        end
        
        try
            Options.clust_method = 'FGK';
            M1a_PrepareDatasets(Options)
        catch
            clearvars -except a year c Options
            close all
            Options.clust_method = 'FCM';
            M1a_PrepareDatasets(Options)
        end
        close all
        if strcmp(Options.clust_method,'FGK')
           clearvars -except a year c Options
           Options.clust_method = 'FCM';
           M1a_PrepareDatasets(Options)
        end
        close all
        clearvars -except a year c Options  
    end
end
fprintf('All patients done - Fuzzy \n')

%% LR models - All

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = true;
Options.Fuzzy = false;
Options.SVM = false;
Options.Kappa = false;
Options.Patients = 'All';

% CV fold
Options.k = 5;

% Partitions (runs)
Options.n_partitions = 10;

% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

for a=1:2
    year = a;
    Options.cluster = [];
    Options.Regression = true;
    Options.Classification = false;
        
    switch year
        case 1
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
             load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...   
                'tbl_in_1yr','var_names1')
            Options.input = tbl_in_1yr;
            Options.year = year;
            Options.var_names = var_names1;
        case 2
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
            load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                'tbl_in_2yr','var_names2')
            Options.input = tbl_in_2yr;
            Options.year = year;
            Options.var_names = var_names2;
    end   
    
    Options.clust_method = [];
    M1a_PrepareDatasets(Options)
 
    close all
    clearvars -except a year c Options  
end

fprintf('All patients done - LR \n')

%% SVM models - All

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = false;
Options.Fuzzy = false;
Options.SVM = true;
Options.Kappa = false;
Options.Patients = 'All';

% CV fold
Options.k = 5;

% Partitions (runs)
Options.n_partitions = 10;

% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

for a=1:2
    year = a;
    Options.cluster = [];
    Options.Regression = true;
    Options.Classification = false;
        
    switch year
        case 1
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
             load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...   
                'tbl_in_1yr','var_names1')
            Options.input = tbl_in_1yr;
            Options.year = year;
            Options.var_names = var_names1;
        case 2
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
            load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                'tbl_in_2yr','var_names2')
            Options.input = tbl_in_2yr;
            Options.year = year;
            Options.var_names = var_names2;
    end   
    
    Options.clust_method = [];
    M1a_PrepareDatasets(Options)
 
    close all
    clearvars -except a year c Options  
end

fprintf('All patients done - SVM \n')

%% Fuzzy Regression models - Metabolic

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = false;
Options.Fuzzy = true;
Options.SVM = false;
Options.Kappa = false;
Options.Patients = 'Metab';

% CV fold
Options.k = 5;

% Partitions (runs)
Options.n_partitions = 10;

% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

for a=1:2
    year = a;
    for c=2:6
        Options.cluster = c;
        Options.Regression = true;
        Options.Classification = false;
        
        switch year
            case 1
                %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
                load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                    'tbl_1yr_metabolic')
                Options.input = tbl_1yr_metabolic;
                Options.year = year;
                Options.var_names = tbl_1yr_metabolic.Properties.VariableNames;
            case 2
                %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
                load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                    'tbl_2yr_metabolic')
                Options.input = tbl_2yr_metabolic;
                Options.year = year;
                Options.var_names = tbl_2yr_metabolic.Properties.VariableNames;
        end
        
        try
            Options.clust_method = 'FGK';
            M1a_PrepareDatasets(Options)
        catch
            clearvars -except a year c Options
            close all
            Options.clust_method = 'FCM';
            M1a_PrepareDatasets(Options)
        end
        close all
        if strcmp(Options.clust_method,'FGK')
           clearvars -except a year c Options
           Options.clust_method = 'FCM';
           M1a_PrepareDatasets(Options)
        end
        close all
        clearvars -except a year c Options  
    end
end
fprintf('Metabolic patients done - Fuzzy \n')

%% LR models - Metabolic

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = true;
Options.Fuzzy = false;
Options.SVM = false;
Options.Kappa = false;
Options.Patients = 'Metab';

% CV fold
Options.k = 5;

% Partitions (runs)
Options.n_partitions = 10;

% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

for a=1:2
    year = a;
    Options.cluster = [];
    Options.Regression = true;
    Options.Classification = false;
        
    switch year
        case 1
            % Load data
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
             load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                  'tbl_1yr_metabolic')
            Options.input = tbl_1yr_metabolic;
            Options.year = year;
            Options.var_names = tbl_1yr_metabolic.Properties.VariableNames;
        case 2
            % Load data
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
            load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                'tbl_2yr_metabolic')
            Options.input = tbl_2yr_metabolic;
            Options.year = year;
            Options.var_names = tbl_2yr_metabolic.Properties.VariableNames;
    end
    
    Options.clust_method = [];
    M1a_PrepareDatasets(Options)
 
    close all
    clearvars -except a year c Options  
end

fprintf('Metabolic patients done - LR \n')

%% SVM models - Metabolic

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = false;
Options.Fuzzy = false;
Options.SVM = true;
Options.Kappa = false;
Options.Patients = 'Metab';

% CV fold
Options.k = 5;

% Partitions (runs)
Options.n_partitions = 10;

% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

for a=1:2
    year = a;
    Options.cluster = [];
    Options.Regression = true;
    Options.Classification = false;
        
    switch year
        case 1
            % Load data
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
             load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                  'tbl_1yr_metabolic')
            Options.input = tbl_1yr_metabolic;
            Options.year = year;
            Options.var_names = tbl_1yr_metabolic.Properties.VariableNames;
        case 2
            % Load data
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
            load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                'tbl_2yr_metabolic')
            Options.input = tbl_2yr_metabolic;
            Options.year = year;
            Options.var_names = tbl_2yr_metabolic.Properties.VariableNames;
    end
    
    Options.clust_method = [];
    M1a_PrepareDatasets(Options)
 
    close all
    clearvars -except a year c Options  
end

fprintf('Metabolic patients done - SVM \n')

%% Fuzzy Regression models - Bariatric 

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = false;
Options.Fuzzy = true;
Options.SVM = false;
Options.Kappa = false;
Options.Patients = 'Bar';

% CV fold
Options.k = 5;

% Partitions (runs)
Options.n_partitions = 10;

% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

for a=1:2
    year = a;
    for c=2:6
        Options.cluster = c;
        Options.Regression = true;
        Options.Classification = false;
        
        switch year
            case 1
                %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
                load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                    'tbl_1yr_bariatric')
                Options.input = tbl_1yr_bariatric;
                Options.year = year;
                Options.var_names = tbl_1yr_bariatric.Properties.VariableNames;
            case 2
                %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
                load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                    'tbl_2yr_bariatric')
                Options.input = tbl_2yr_bariatric;
                Options.year = year;
                Options.var_names = tbl_2yr_bariatric.Properties.VariableNames;
        end
        
        try
            Options.clust_method = 'FGK';
            M1a_PrepareDatasets(Options)
        catch
            clearvars -except a year c Options
            close all
            Options.clust_method = 'FCM';
            M1a_PrepareDatasets(Options)
        end
        close all
        if strcmp(Options.clust_method,'FGK')
           clearvars -except a year c Options
           Options.clust_method = 'FCM';
           M1a_PrepareDatasets(Options)
        end
        close all
        clearvars -except a year c Options  
    end
end
fprintf('Bariatric patients done - Fuzzy \n')

%% LR models - Bariatric

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = true;
Options.Fuzzy = false;
Options.SVM = false;
Options.Kappa = false;
Options.Patients = 'Bar';

% CV fold
Options.k = 5;

% Partitions (runs)
Options.n_partitions = 10;

% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

for a=1:2
    year = a;
    Options.cluster = [];
    Options.Regression = true;
    Options.Classification = false;
        
    switch year
        case 1
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
             load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                'tbl_1yr_bariatric')
            Options.input = tbl_1yr_bariatric;
            Options.year = year;
            Options.var_names = tbl_1yr_bariatric.Properties.VariableNames;
        case 2
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
            load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                'tbl_2yr_bariatric')
            Options.input = tbl_2yr_bariatric;
            Options.year = year;
            Options.var_names = tbl_2yr_bariatric.Properties.VariableNames;
    end
    
    Options.clust_method = [];
    M1a_PrepareDatasets(Options)
 
    close all
    clearvars -except a year c Options  
end

fprintf('Bariatric patients done - LR \n')

%% SVM models - Bariatric

% Select Model options
Options.Polynomial = false; % Default 1st order, change in SETUP MODEL ALGORITHM section
Options.Logistic = false;
Options.Fuzzy = false;
Options.SVM = true;
Options.Kappa = false;
Options.Patients = 'Bar';

% CV fold
Options.k = 5;

% Partitions (runs)
Options.n_partitions = 10;

% Select just one metric
Options.TWL = true;     % Total weight loss 'TWL'
Options.EBWL = false;   % Excess body weight loss 'EBWL'
Options.EBMIL = false;  % Excess BMI loss 'EBMIL'
Options.WL = false;     % Absolute weight loss 'WL'

for a=1:2
    year = a;
    Options.cluster = [];
    Options.Regression = true;
    Options.Classification = false;
        
    switch year
        case 1
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
             load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                'tbl_1yr_bariatric')
            Options.input = tbl_1yr_bariatric;
            Options.year = year;
            Options.var_names = tbl_1yr_bariatric.Properties.VariableNames;
        case 2
            %load('(2018-07-24 08 57 26) CZE(1yr&2yr)_reg_69var.mat',...
            load('(2019-03-17 01 03 42) CZE(1yr&2yr)_reg_69var.mat',...
                'tbl_2yr_bariatric')
            Options.input = tbl_2yr_bariatric;
            Options.year = year;
            Options.var_names = tbl_2yr_bariatric.Properties.VariableNames;
    end
    
    Options.clust_method = [];
    M1a_PrepareDatasets(Options)
 
    close all
    clearvars -except a year c Options  
end

fprintf('Bariatric patients done - SVM \n')
