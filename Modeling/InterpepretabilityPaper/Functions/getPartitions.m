function [labels,cvp1,data_train,data_test,cvp2] = getPartitions...
    (data_table, filename_partitions, Options)
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