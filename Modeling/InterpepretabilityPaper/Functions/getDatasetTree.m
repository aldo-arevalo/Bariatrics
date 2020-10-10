function [cvp1,data_train,data_test,cvp2,Options] = getDatasetTree(Options)
close all
try
   Options.input.vID120836 = [];
catch
   disp('vID120836 NOT REMOVED')
end

try
    Options.input.hypert = [];
catch
    disp('Hypertension NOT REMOVED')
end

try
    Options.input.dyslip = [];
catch
    disp('dyslipidemia NOT REMOVED')
end

try
    Options.input.Age_ = [];
catch
    disp('Age_ NOT REMOVED')
end

%% Options

if Options.Classification == true
    tmp_thr_test = exist('thr_twl','var');
    % Define threshold for TWL from 0 - 1
    if tmp_thr_test == 0
       thr_twl = 0.25;
    end
    clearvars tmp_thr_test
end

%% Load data

% Columns::
    %   + Partient code
    %   + Inputs: features (2:end-1)
    %   + Output: TWL

% Arranje tables
if Options.year == 1
    var_names = Options.input.Properties.VariableNames;
    data_table = Options.input;
    
    % Others metrics
    if Options.EBWL == true
        load('Classifiers2017+-3months.mat', 'Classifiers1yr')
        data_table.TWL = [];
        data_table_2 = Classifiers1yr;
        data_table_2(:,{'EBMIL','WL','TWL'}) = [];
        data_table = innerjoin(data_table_2,data_table,...
        'LeftKeys','PatientCode','RightKeys','PatientCode');
    end
    if Options.EBMIL == true
       load('Classifiers2017+-3months.mat', 'Classifiers1yr')
       data_table.TWL = [];
       data_table_2 = Classifiers1yr;
       data_table_2(:,{'EBWL','WL','TWL'}) = [];
       data_table = innerjoin(data_table_2,data_table,...
        'LeftKeys','PatientCode','RightKeys','PatientCode');
    end
    if Options.WL == true
       load('Classifiers2017+-3months.mat', 'Classifiers1yr')
       data_table.TWL = [];
       data_table_2 = Classifiers1yr;
       data_table_2(:,{'EBWL','EBMIL','TWL'}) = [];
       data_table = innerjoin(data_table_2,data_table,...
        'LeftKeys','PatientCode','RightKeys','PatientCode');
    end
else % year == 2
    var_names = Options.input.Properties.VariableNames;
    data_table = Options.input;

    if Options.EBWL == true
       load('Classifiers2017+-3months.mat', 'Classifiers2yr')
       data_table.TWL = [];
       data_table_2 = Classifiers2yr;
       data_table_2(:,{'EBMIL','WL','TWL'}) = [];
       data_table = innerjoin(data_table_2,data_table,...
            'LeftKeys','PatientCode','RightKeys','PatientCode');
    end
    if Options.EBMIL == true
       load('Classifiers2017+-3months.mat', 'Classifiers2yr')
       data_table.TWL = [];
       data_table_2 = Classifiers2yr;
       data_table_2(:,{'EBWL','WL','TWL'}) = [];
       data_table = innerjoin(data_table_2,data_table,...
            'LeftKeys','PatientCode','RightKeys','PatientCode');
    end
    if Options.WL == true
        load('Classifiers2017+-3months.mat', 'Classifiers2yr')
       data_table.TWL = [];
       data_table_2 = Classifiers2yr;
       data_table_2(:,{'EBWL','EBMIL','TWL'}) = [];
       data_table = innerjoin(data_table_2,data_table,...
            'LeftKeys','PatientCode','RightKeys','PatientCode');
    end
end

idx_missing = any(ismissing(data_table),2);
data_table(idx_missing,:) = [];

% Select output and input
if Options.TWL == true
    metric='TWL';
end
if Options.EBWL == true
  var_names(strcmp(var_names,'TWL')) = [];
%   Y = data_table.EBWL;
  data_table = data_table(:,var_names);
  metric='EWL';
end
if Options.EBMIL == true
  var_names(strcmp(var_names,'TWL')) = [];
%   Y = data_table.EBMIL;  
  data_table = data_table(:,var_names);
  metric='EBMIL';
end
if Options.WL == true
  var_names(strcmp(var_names,'TWL')) = [];
%   Y = data_table.WL;  
  data_table = data_table(:,var_names);
  metric='WL';
end

Options.metric = metric;
data_table.PatientCode = [];
   
if Options.Regression == true
    fprintf('Regression')
else % Classification == true
    total_sucess = sum(data_table.(metric));
    total_non_sucess = size(data_table.(metric),1)- total_sucess;
    Options.Class1 = total_sucess;
    Options.Class0 = total_non_sucess;
end

%% Convert numerical vars to categorical

% Did patient is treated by a specialist during last year? - Hist Dr
tmp_exist_vID120883 = strcmp(data_table.Properties.VariableNames,...
    'vID120883');
if sum(tmp_exist_vID120883) == 1
    data_table = convertvars(data_table,{'vID120883'},'categorical');
end

% Did patient has been treated by a doctor for high blood pressure? - Hist Hyper
tmp_exist_vID120848 = strcmp(data_table.Properties.VariableNames,...
    'vID120848');
if sum(tmp_exist_vID120848) == 1
    data_table.vID120848 = data_table.vID120848 > 0.2;
    data_table = convertvars(data_table,{'vID120848'},'categorical');
end

% Did patient come due to other physical complaints? - Hist Other
tmp_exist_vID120881 = strcmp(data_table.Properties.VariableNames,...
    'vID120881');
if sum(tmp_exist_vID120881) == 1
    data_table.vID120881 = data_table.vID120881 > 0.2;
    data_table = convertvars(data_table,{'vID120881'},'categorical');
end

% Gender, roken, diabet, hypert, dyslip
data_table = convertvars(data_table,{'Sex','roken','diabet'},...
    'categorical');

tmp_exist_Sex = strcmp(data_table.Properties.VariableNames,'Sex');
tmp_exist_roken = strcmp(data_table.Properties.VariableNames,'roken');
tmp_exist_diabet = strcmp(data_table.Properties.VariableNames,'diabet');

% Identify categorical variables
tmp_idx_categorical = [tmp_exist_vID120883; ...
    tmp_exist_vID120848; tmp_exist_vID120881; tmp_exist_Sex; ...
    tmp_exist_roken; tmp_exist_diabet];
    
Options.isCategoricalPredictor = logical(sum(tmp_idx_categorical));

clearvars tmp_exist_vID120881 tmp_exist_vID120848 tmp_exist_vID120883 ...
    tmp_idx_categorical tmp_exist_Sex tmp_exist_roken tmp_exist_diabet
%% Generate Hold-out
if Options.Polynomial == 0
    
    % Create folder Paritions in case it has not been created
    tmp_folder_exist = exist('Partitions','file');
    
    if ne(tmp_folder_exist, 7)
        mkdir 'Partitions'
    end
    
    % Create file for naming the files that contains the partitions
    if Options.Regression == true
        filename_partitions = (['Partition_',Options.JOB,'_Year',...
            char(string(Options.year)),'_Reg_',metric,'_k',...
            char(string(Options.k)),Options.Patients,'.mat']);
    else % Classification
        filename_partitions = (['Partition_',Options.JOB,'_Year',...
            char(string(Options.year)),'_',char(string(thr_twl*100)),...
            metric,'_k',char(string(Options.k)),'_',...
            Options.Patients,'.mat']);
    end
    
    tmp_file_exist = exist(['Partitions/',filename_partitions],'file');

    if ne(tmp_file_exist,2) == true
       [labels,cvp1,data_train,data_test,cvp2] = getPartitions...
           (data_table, filename_partitions, Options);
    else
        load(['Partitions/',filename_partitions],...
            'labels','cvp1','data_train','data_test','cvp2')
    end
    
    Options.Paritions_file = filename_partitions;
    
    % Change variable names for plotting
    Options.var_names_modif = labels;
    
    tmp_idx_vars = strcmp(data_table.Properties.VariableNames, ...
        Options.metric);
    Options.predictorNames = labels;
    Options.isCategoricalPredictor = Options.isCategoricalPredictor(~tmp_idx_vars);
    data_table.Properties.VariableNames(~tmp_idx_vars) = labels;
    data_train.Properties.VariableNames(~tmp_idx_vars) = labels;
    data_test.Properties.VariableNames(~tmp_idx_vars) = labels;
    Options.data_table = data_table;
    
end