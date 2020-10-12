function [holdPartition,trainingData,testData,...
    cvPartition,Options] =  getInputSubsets(Options, filename)

% Function to obtain the cross validation partition for training/validation
% and hold out as defined in 'Options' structure.
%
% INPUTS:
%       Options:  Structure that contains the general settings.
%       filename: String that contains the path to the curated dataset.
%       year:     Integer specifying the follow-up year to predict.
%
% OUTPUTS:
%       holdPartition: Stratified hold out partition.
%       trainingData:  Stratified cross validation partition.
%       trainingData:  Data set for training (x features, y output variable).
%       testData:      Data for testing (x features, y output variable).
%
% Dependencies: M1_Modeling
%
% Author: Aldo Arévalo
% Date: 10/11/2020

switch Options.year
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
[holdPartition,trainingData,testData,cvPartition,Options] = getDataset(Options);