%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dataset: Catharina Hospital
% Language: MATLAB 2016a and forward
% Author: Aldo R. Arevalo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run scripts in correct order and save variables
clear all

% Create a folder where tables are going to be saved
d=datestr(date);
subfold = 'CZE_data_extraction';

A = exist(subfold,'dir');

if A == 0
    mkdir(subfold,d)
    mkdir([subfold,'/',d,'/'],'Plots')
end

clearvars A

%% Import Excel spreadsheets into MATLAB
    run S1_ImportData2018
    close all
    save([subfold,'/',d,'/ImportedTables']);
    fprintf('Stage 1 completed\n')
    
%% Z-transform for Lab analyses variables and Age
    SD = 2;    
    run S2a_RefValuesLab
    close all
    save([subfold,'/',d,'/LabZt']);
    fprintf('Stage 2a completed\n')   

%% Create a Structure for Lab Results
    run S2b_CalcLabResultsPerPatient2018
    close all
    save([subfold,'/',d,'/LabResultsStructure']);
    fprintf('Stage 2b completed\n')  

%% Create a Structure for Physical Measures, includes QoL info
    
    % Minimal frequency for variables
    a=10; % Fuchs, C. used 250
    run S3_CalcPhMeasPerPatient2018
    close all
    save([subfold,'/',d,'/PhMsStructure']);
    clearvars b ColumSruct Creatinine_meas LabResultsPerPatient myStruct ...
        tbl2 analysisList_empty
    fprintf('Stage 3 completed\n')

%% Create a Structure for Pre-Surgical Screening
    run S4_CalcPPOSPerPatient2018
    close all
    clearvars A b B vrID PPOSStruct PPOSPerPatient
    save([subfold,'/',d,'/PPOSstructure']);
    fprintf('Stage 4 completed\n')

%% Create tables with minima, maxima, last measured, mean and number of times measured: Lab Results
    run S5_TemporarySeriesLab2018
    close all
    save([subfold,'/',d,'/LabResultsTimeSeries']);
    fprintf('Stage 5 completed\n')

%% Create tables with minima, maxima, last measured, mean and number of
    % timesmeasured: Physical Measures
    SD = 2;
    run S6_TemporarySeriesPhMs2018
    save([subfold,'/',d,'/PhMsTimeSeries']);
    fprintf('Stage 6 completed\n')
    clearvars -except d subfold

%% Create tables with minima, maxima, last measured, mean and number of
    % timesmeasured: Pre-Surgical Screening
    run S7_TemporarySeriesPPOS2018
    save([subfold,'/',d,'/PPOSTimeSeries']);
    fprintf('Stage 7 completed\n')
    clearvars -except d subfold

%% Extract pre-operative variables: weight, height, BMI, BRI, ABSI, TBFM 
    run S8_PreOperativeValues2018
    save([subfold,'/',d,'/PreOperativeValues']);
    fprintf('Stage 8 completed\n')
    clearvars -except d subfold

%% Estimate Classes 
    run S9_Classifier2018
    save([subfold,'/',d,'/Classifiers2017+-3months']);
    fprintf('Stage 9 completed\n')
    clearvars -except d subfold
    
%% Build input table
% Select as desired. The default options is for building regression models.

% Missing data. Mantain missing data?
% True = Mantain missing data (may need an imputation method)
% False = List-wise deletion method
missing_data = false;
    
if missing_data == false % Apply deletion method in input matrix
    % Define bias
        % Deleting all variables below threshold = 0
        % Deleting all variables except WC, BRI and ABSI = 1
        bias = 0;
end

% Binary classification?
bin_class = false;
    if bin_class == true
        % Define threshold for TWL from 0 - 1
        thr_twl = 0.25;
        id_twl = [num2str(thr_twl*100),'%TWL'];      
    end
% Generate histograms of distribution of success class among gender,
% type of surgery, BMI and comorbidities
hist_dist = false;

% Multivariate classification?
multi_class = false;
    if multi_class == true
        % Class 1: 0-25%TWL
        class0 = 0.25;
        % Class 2: 25-35%TWL
        class1 = 0.35;
        % Class 3: 35-45%TWL
        class2 = 0.45;
        % Class 4: >45%TWL
    end

% Regression model?
reg_model = true;

% Delete all binary features?
bin_delete = false;

% Create Quality of Life score?
QoL = true;
 % If true, all these binary varaibles will be deleted.

% Save figures or plots?
keep_fig = false;
    
% Create new folder?
create_folder = true;

% Run code
dd = d;
run S10_BuildingInputTable2018
d = dd;
save([subfold,'/',d,'/FinalTables']);
close all

%% Finishing up
fprintf('Extraction of data finished \n')
beep
