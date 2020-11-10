% Run scripts in correct order and save variables
clearvars

% Create folder where the data is going to be stored
d=datestr(date);
subfold = 'CZE_data_extraction';
A = exist([subfold,'/',d],'dir');

% Folder for storing plots and for raw matrixes (mat files)
if A == 0
    mkdir(subfold,d)
    mkdir([subfold,'/',d,'/'],'Plots')
    mkdir([subfold,'/',d,'/'],'RawData')
    addpath(genpath([subfold,'/',d,'/RawData']))
end
clearvars A

% Folder where raw data (xlsx files) are stored
addpath('RawData')

% Scripts called by MainScript
addpath('SecondaryScripts')

% Functions needed by the secondary scripts
addpath('Functions')
%% Import Excel spreadsheets into MATLAB
    run S1_ImportData
    close all
    save([subfold,'/',d,'/ImportedTables'],'FollowUpData','GeneralData',...
        'LabResults','PhMsData','PPOSdata','ScreeningData',...
        'SurgeryData','total_patientsID','ver','d','subfold');
    fprintf('Stage 1 completed\n')
    
    clearvars -except d subfold
%% Z-transform for Lab analyses variables and Age
    SD = 2;    
    run S2a_RefValuesLab
    close all
    save([subfold,'/',d,'/LabZt']);
    fprintf('Stage 2a completed\n')   

    clearvars -except d subfold
%% Create a Structure for Lab Results
    run S2b_CalcLabResultsPerPatient
    close all
    save([subfold,'/',d,'/LabResultsStructure']);
    fprintf('Stage 2b completed\n')  
    
    clearvars -except d subfold
%% Create a Structure for Physical Measures, includes QoL info
    
    % Minimal frequency for variables
    a=10; % Fuchs, C. used 250
    run S3_CalcPhMeasPerPatient
    close all
    save([subfold,'/',d,'/PhMsStructure'],'d','PhMsData',...
        'PhMsPerPatient','subfold','tbl2');
    clearvars b ColumSruct Creatinine_meas LabResultsPerPatient myStruct ...
        tbl2 analysisList_empty
    fprintf('Stage 3 completed\n')

    clearvars -except d subfold a
%% Create a Structure for Pre-Surgical Screening
    run S4_CalcPPOSPerPatient
    close all
    clearvars A b B vrID PPOSStruct PPOSPerPatient
    save([subfold,'/',d,'/PPOSstructure']);
    fprintf('Stage 4 completed\n')
    
    clearvars -except d subfold
%% Create tables with minima, maxima, last measured, mean and number of times measured: Lab Results
    run S5_TemporarySeriesLab
    close all
    save([subfold,'/',d,'/LabResultsTimeSeries'],'LabLastPatient',...
        'LabMaxPatient','LabMeanPatient','LabMedianPatient',...
        'LabMinPatient','LabNMPatient','LabNTPrior',...
        'LabPackTimesPatient','LabTotalNT','total_patientsID',...
        'var_list','lab_pack');
    fprintf('Stage 5 completed\n')
    
    clearvars -except d subfold
%% Create tables with minima, maxima, last measured, mean and number of
    % timesmeasured: Physical Measures
    SD = 2;
    run S6_TemporarySeriesPhMs
    close all
    save([subfold,'/',d,'/PhMsTimeSeries'],'PhMsData','PhMsLastPatient',...
        'PhMsMedianPatient','PhMsNMPatient','PhNTPrior','PhTotalNT',...
        'subfold','d');
    fprintf('Stage 6 completed\n')
    
    clearvars -except d subfold
%% Create tables with minima, maxima, last measured, mean and number of
    % timesmeasured: Pre-Surgical Screening
    run S7_TemporarySeriesPPOS
    close all
    save([subfold,'/',d,'/PPOSTimeSeries'],'d','subfold','tbl3',...
        'PPOSLastPatient','PPOSMedianPatient');
    fprintf('Stage 7 completed\n')
    
    clearvars -except d subfold
%% Extract pre-operative variables: weight, height, BMI, BRI, ABSI, TBFM 
    run S8_PreOperativeValues
    close all
    save([subfold,'/',d,'/PreOperativeValues'],'PreOpLast','PreOpMax',...
        'PreOpMean','PreOpMedian','PreOpMin','PreOpNMTotal',...
        'patients_list','bmi_max','weight_max','height_mean');
    fprintf('Stage 8 completed\n')
    
    clearvars -except d subfold
%% Estimate Classes 
    run S9_Classifier
    close all
    save([subfold,'/',d,'/Classifiers',d(8:11),'+-3months'],'subfold',...
        'd','Classifiers1yr','Classifiers2yr');
    fprintf('Stage 9 completed\n')
    
    clearvars -except d subfold
%% Build input table

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
bin_class = true;
    if bin_class == true
        % Define threshold for TWL from 0 - 1
        thr_twl = 0.25;
        id_twl = [num2str(thr_twl*100),'%TWL'];      
    end
% Generate histograms of distribution of success class among gender,
% type of surgery, BMI and comorbidities
hist_dist = true;

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
reg_model = false;

% Delete all binary features?
bin_delete = false;

% Create Quality of Life score?
QoL = true;
 % If true, all these binary varaibles will be deleted.

% Save figures?
keep_fig = true;
    
% Create new folder?
create_folder = true;

% Run code
dd = d;
run S10_BuildingInputTable

close all

fprintf('Stage 10 completed\n')

clearvars -except dd subfold filename
%% PreOperative values for Raw data
d = dd;

run S11_PreOperativeValuesRAW
close all
save([subfold,'/',d,'/PreOperativeValuesRAW'],'bmi_max','BMIall',...
    'd','DBMI','DHeight','DWC','DWeight','HeightALL',...
    'LabPackTimesPatient','LabResultsPerPatient','list_PatientID',...
    'patients_list','stat','WeightALL','weight_max','mWeight','sWeight',...
    'WCall');
fprintf('Stage 11 completed\n')

clearvars -except d subfold filename
    
%% Estimate TWL for RAW data

run S12_ClassifierRAW
fprintf('Stage 12 completed\n')

clearvars -except d subfold filename

%% Plotting raw data

run S13_Rawdata_analysis
close all
save([subfold,'/',d,'/RawPlots&ContingencyTables'],'tbl_cross1yr',...
    'tbl_cross2yr','tbl_h_1yr','tbl_h_2yr','tbl_pval_1yr','tbl_pval_2yr',...
    'RAW_1yr','RAW_2yr');
fprintf('Stage 13 completed\n')

%% Finishing up
fprintf('Extraction of data finished \n')
beep