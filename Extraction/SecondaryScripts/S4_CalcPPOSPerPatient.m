%% Find Pre-Operative measurements per patient
% In Dutch: Polikliniek Pre-Operatieve Screening
% In English: outpatient department of pre-surgical screening.
% Make unique list U with measurement codes vraagID which
% have a frequency >=a

A = exist('PPOSdata','var');

if A == 0
    load([subfold,'/',d,'/ImportedTables.mat'],'PPOSdata')
end

% Obtain statistics
vrID = table2cell(PPOSdata(:,'questionID'));
x=unique(PPOSdata.questionID);
tbl3={};

for i=1:length(x)
    var = x(i);
    idx_vrID = strcmp(vrID, var);
    
    tbl3(i,1)= cellstr(var);
    tbl3(i,2)= num2cell(sum(idx_vrID));
end

% tbl3=cell2table(tbl3);

% Specify the height of the table and delete measurementes below 100 times
idx_meas_below_a_times = [tbl3{:,2}]' < a;

U = tbl3(~idx_meas_below_a_times,1);

clear i a x var idx_vrID idx_meas_below_a_times
%% Initialize Structure
% Extract vraagIDs in a str cell array
var_names = unique(U);

data_struct.Result = {};
data_struct.days_from_ok = NaN;

% Create front row with Question ID (vraagID) names
PPOSStruct.PatientCode = NaN;
for i = 1:length(var_names)
    PPOSStruct.(var_names{i}) = data_struct;
end

clearvars i U var_names
%% Fill in the structure for each patient

list_PatientID = unique(PPOSdata.PatientCode);

% Pre-allocate structs for all patients
PPOSPerPatient = repmat(PPOSStruct,length(list_PatientID),1);

for i=1:length(list_PatientID)
    patient_ID = list_PatientID(i);
    
    PPOSPerPatient(i).PatientCode = patient_ID;
    
    idx_patient_ID = PPOSdata.PatientCode == patient_ID;
    
    s = PPOSdata(idx_patient_ID,:);
    
    var_names = unique(s.questionID);
    
    for j = 1:length(var_names)
        var_name = var_names(j);
        idx_var = strcmp(s.questionID, var_name);
        
        vars = table2cell(s(idx_var,5:7));
        n = size(vars,1);
        
        variables_list = {};
        for k = 1:n
            % Check for strings and NaNs
            idx_is_not_char = find(~cellfun(@ischar,vars(k,:)));
            idx_valid_meas = idx_is_not_char(~cellfun(@isnan, vars(k, idx_is_not_char)));
            
            variables_list{k} = cell2mat(vars(k,idx_valid_meas));
        end
        days_from_ok_list = s.days_from_ok(idx_var);
        
        PPOSPerPatient(i).(var_name{1,1}).Result = variables_list;
        PPOSPerPatient(i).(var_name{1,1}).days_from_ok = days_from_ok_list;
    end
    clear idx_var vars n s var_names idx_patient_ID var_name ...
        days_from_ok_list k j i idx_is_not_char idx_valid_meas patient_ID
end

% The question ID 120864 will be removed
% Exclusion done in PstgreSQL
% PPOSPerPatient(:,'120864') = [];
%% Histogram
% Variables or question IDs extracted
C = sortrows(tbl3,2);
figure('Name','Questions frequency')
b=bar(cell2mat(C(:,end)));
xlabel('Question ID');
ylabel('Frequency');
print([subfold,'/',d,'/Plots/PPOSs4'],'-deps')
print([subfold,'/',d,'/Plots/PPOSs4c'],'-depsc')
print([subfold,'/',d,'/Plots/PPOSs4'], '-dpdf', '-fillpage')

clear data_struct variables_list C tbl3
