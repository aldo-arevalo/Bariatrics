function [labels,cvp1,data_train,data_test,cvp2] = getPartitions...
    (data_table, filename_partitions, Options)
% Function that creates the training/validation and testing subsets. Also
% creates the labels for comprehensive reading of plots.
%
% INPUTS:
%       data_table:          Input dataset (table).
%       filename_partitions: String containing the path to save files.
%       Options:             Structure that contains the general settings.
%
% OUTPUTS:
%       labels:     Vector containing the feature names for printing.
%       cvp1:       Stratified hold out partition.
%       data_train: Data set for training (x features, y output variable).
%       data_test:  Data for testing (x features, y output variable).
%       cvp2:       Stratified cross validation partition.
%
% Dependencies: getDataset
%
% Author: Aldo Arévalo
% Date: 10/11/2020

% Generate labels for plots
labels = getVarsLabels(data_table.Properties.VariableNames,Options);

%% Make sperate training and test sets

cvp1 = cvpartition(data_table.(Options.metric),'HoldOut',Options.Hold_out,...
    'Stratify',true);

data_train = data_table(cvp1.training,:);

data_test = data_table(cvp1.test,:);

%% CV partition and save

cvp2 = cvpartition(data_train.(Options.metric), 'kfold', Options.k, ...
    'Stratify', true);

save(['Partitions/',filename_partitions],...
    'labels', 'cvp1', 'data_train', 'data_test', 'cvp2')