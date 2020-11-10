%% Build Table for time series variables
% Mean, median, maximum, minimum, last value before surgery

% Load data
A = exist('GeneralData','var');
if A == 0
    load([subfold,'/',d,'/ImportedTables.mat'],'GeneralData')
end

A = exist('PhMsData','var');
if A == 0
    load([subfold,'/',d,'/ImportedTables.mat'],'PhMsData')
end

A = exist('LabPackTimesPatient','var');
if A == 0
    load([subfold,'/',d,'/LabResultsTimeSeries.mat'], 'LabPackTimesPatient')
end

clearvars A
%% Delete BMI, HEIGHT, WAIST CIRCUMF. and WEIGHT variables
% These variables will be used to estimate pre-operative features
    % Weight (gewicht)
    PhMsData(ismissing(PhMsData.vraagID, 111972),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 112208),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 113882),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 117529),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 118486),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 121008),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 121101),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 124719),:)=[];
    
    % BMI
    PhMsData(ismissing(PhMsData.vraagID, 116086),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 118678),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 120703),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 124718),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 118488),:)=[];
    
    % Height (Lengte)
    PhMsData(ismissing(PhMsData.vraagID, 116084),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 118487),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 118677),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 121100),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 132959),:)=[];
    
    % Waist Circumference (buikomvang)
    PhMsData(ismissing(PhMsData.vraagID, 118489),:)=[];
    
%% Gather variables names

% Columns/variables: Convert vraagID into elegible feature names
questionID = cell(size(PhMsData,1),1);

for i=1:size(PhMsData,1)
    questionID{i} = ['qID' num2str(PhMsData.vraagID(i))];
end

% Columns/variables: Concatenate questionID colum to PhMsData
PhMsData.questionID=questionID;

% Obtain statistics from data gathered before surgery
tbl2=tabulate(PhMsData.questionID);

% Minimal frequency
a=10; % Fuchs, C. used 250

% Specify the height of the table and delete measurementes below 100 times
idx_meas_below_a_times = [tbl2{:,2}]' < a;

U = tbl2(~idx_meas_below_a_times,1);

clearvars a questionID

%% Build tables

% Extract vraagID
qIDs = unique(U);

% Rows
u = unique(GeneralData.PatientNr);

% Normal blood pressure in adults is lower than 
% 120/80 mmHg. Hypotension is blood pressure that's
% lower than 90/60 mmHg. National Heart, Lung, and
% Blood Institute. Z-transform
sSys = (120 - 90)/4;     % sigma = (RefUp - RefLow)/4
mSys = 90 + SD*sSys;     % mu = RefLow + 2s

sDias = (80 - 60)/4;     % sigma = (RefUp - RefLow)/4
mDias = 60 + SD*sDias;   % mu = RefLow + 2s

% Pre-allocate matrix
prealloc_num = array2table(NaN(length(u), length(qIDs)+3));
prealloc_cell = array2table(repmat({NaN},length(u), length(qIDs)+3));


% PhMsMeanPatient = prealloc_num;
%PhMsMinPatient = prealloc_num;
%PhMsMaxPatient = prealloc_num;
PhMsNMPatient = prealloc_num;
PhMsLastPatient = prealloc_cell;
PhMsMedianPatient = prealloc_num;
PhNTPrior = prealloc_num;
PhTotalNT = prealloc_num;

% Add column names
var_list = vertcat('PatientCode',qIDs,'Systolic','Diastolic');
% PhMsMeanPatient.Properties.VariableNames = var_list;
%PhMsMinPatient.Properties.VariableNames = var_list;
%PhMsMaxPatient.Properties.VariableNames = var_list;
PhMsNMPatient.Properties.VariableNames = var_list;
PhMsLastPatient.Properties.VariableNames = var_list;
PhMsMedianPatient.Properties.VariableNames = var_list;
PhNTPrior.Properties.VariableNames = var_list;
PhTotalNT.Properties.VariableNames = var_list;

for i=1:length(u)
    patient_ID = u(i);
    idx_patient_ID = [PhMsData.PatientCode] == patient_ID;
    id_patient_screen = table2array(LabPackTimesPatient(:,{'PatientCode'})) == patient_ID;
    screen_day = table2array(LabPackTimesPatient(id_patient_screen,{'ScreenDay'}));
    max_screen_day = screen_day - 15;
    min_screen_day = screen_day + 15;
    
%     PhMsMeanPatient.PatientCode(i) = patient_ID;
    %PhMsMinPatient.PatientCode(i) = patient_ID;
    %PhMsMaxPatient.PatientCode(i) = patient_ID;
    PhMsNMPatient.PatientCode(i) = patient_ID; % During screen time window
    PhMsLastPatient.PatientCode(i) = {patient_ID};
    PhMsMedianPatient.PatientCode(i) = patient_ID;
    PhNTPrior.PatientCode(i) = patient_ID; % Before surgery screen up to 1 year
    PhTotalNT.PatientCode(i) = patient_ID; % Before 1 year
    
    s = PhMsData(idx_patient_ID,:);
    in_windowX = arrayfun(@(x) x >= -365, s.days_from_ok);
    
    for j = 1:length(qIDs)
        var_name = qIDs{j};
        idx_var = strcmp(s.questionID,var_name);
        
        if sum(idx_var)>= 1
            in_twindow = arrayfun(@(x) x >= -365, table2array(s(idx_var,{'days_from_ok'})));
            
            if sum(in_twindow) >= 1
                % Combine indexes
                idx_conj = sum(horzcat(idx_var,in_windowX),2)==2;
                
                Ant = table2cell(s(idx_conj,{'antwoord'}));
                Ant = cellfun(@str2double,Ant);
                varstmp = table2array(s(idx_conj,{'value01','value02','days_from_ok'}));
                vars = horzcat(Ant,varstmp);
                
                % vars(isnan(vars)) = [];
                n = size(vars,1);

                variables_list = {};
                for k = 1:n
                    tmp = vars(k,:);
                    tmp(isnan(tmp)) = [];
                    variables_list{k} = tmp;
                end
                clearvars k tmp Ant varstmp
                    
                if isempty(variables_list{1,1}) == 0
                    % Id of tension; test patient 14629437870
                    if strcmp(var_name,'qID119058') == 1
                       % Find last measurment before surgery
                       vars(:,sum(isnan(vars))>=size(vars,1)) = [];
                       var_time = sortrows(vars, size(vars,2));
                       row_nans = any(~isnan(var_time(:,1:2)),2);
                        
                       % Systolic and Diastolic                      
                        if size(var_time,1) == 1
                           var_time(isnan(var_time)) = [];

                           zSys = (max(var_time(1,1:2))-mSys)/sSys;
                           
                           PhMsLastPatient.Systolic(i) = num2cell(zSys);
                           PhMsMedianPatient.Systolic(i) = zSys;
                           %PhMsMaxPatient.Systolic(i) = zSys;
                           %PhMsMinPatient.Systolic(i) = zSys;
%                            PhMsMeanPatient.Systolic(i) = zSys;
                           
                           zDias = (min(var_time(1,1:2))-mDias)/sDias;
                           
                           PhMsLastPatient.Diastolic(i) = num2cell(zDias);
                           PhMsMedianPatient.Diastolic(i) = zDias;
                           %PhMsMaxPatient.Diastolic(i) = zDias;
                           %PhMsMinPatient.Diastolic(i) = zDias;
%                            PhMsMeanPatient.Diastolic(i) = zDias;
                           
                           % Check if measurement was taken during screening
                           PhMsNMPatient.Systolic(i) = sum(arrayfun(@(x) x >= max_screen_day & x <= min_screen_day , var_time(:,3)));
                           PhMsNMPatient.Diastolic(i) = sum(arrayfun(@(x) x >= max_screen_day & x <= min_screen_day , var_time(:,3)));
                           % Chech if measurement was taken before screening
                           % but up to a year
                           PhNTPrior.Systolic(i) = sum(arrayfun(@(x) x >= -365 & x <= max_screen_day , var_time(:,3)));
                           PhNTPrior.Diastolic(i) = sum(arrayfun(@(x) x >= -365 & x <= max_screen_day , var_time(:,3)));
                           % Check if measurement was taken before a year of
                           % surgery
                           PhTotalNT.Systolic(i) = length(var_time(:,3));
                           PhTotalNT.Diastolic(i) = length(var_time(:,3));
                        else % More rows, more observations
                           syst = NaN(size(var_time,1),1);
                           dias = NaN(size(var_time,1),1);

                           for p=1:size(var_time,1)
                               syst(p,1) = (max(var_time(p,1:2))-mSys)/sSys;
                               dias(p,1) = (min(var_time(p,1:2))-mDias)/sDias;
                           end
                           clearvars p

                           syst = syst(row_nans,1);
                           dias = dias(row_nans,1);

                           PhMsLastPatient.Systolic(i) = num2cell(syst(end,1));
                           PhMsMedianPatient.Systolic(i) = median(syst,'omitnan');
                           %PhMsMaxPatient.Systolic(i) = max(syst);
                           %PhMsMinPatient.Systolic(i) = min(syst);
%                            PhMsMeanPatient.Systolic(i) = mean(syst,'omitnan');

                           PhMsLastPatient.Diastolic(i) = num2cell(dias(end,1));
                           PhMsMedianPatient.Diastolic(i) = median(dias,'omitnan');
                           %PhMsMaxPatient.Diastolic(i) = max(dias);
                           %PhMsMinPatient.Diastolic(i) = min(dias);
%                            PhMsMeanPatient.Diastolic(i) = mean(dias,'omitnan');
                           
                           % Check if measurement was taken during screening
                           PhMsNMPatient.Systolic(i) = sum(arrayfun(@(x) x >= max_screen_day & x <= min_screen_day , var_time(row_nans,3)));
                           PhMsNMPatient.Diastolic(i) = sum(arrayfun(@(x) x >= max_screen_day & x <= min_screen_day , var_time(row_nans,3)));
                           % Chech if measurement was taken before screening
                           % but up to a year
                           PhNTPrior.Systolic(i) = sum(arrayfun(@(x) x >= -365 & x <= max_screen_day , var_time(row_nans,3)));
                           PhNTPrior.Diastolic(i) = sum(arrayfun(@(x) x >= -365 & x <= max_screen_day , var_time(row_nans,3)));
                           % Check if measurement was taken before a year of
                           % surgery
                           PhTotalNT.Systolic(i) = length(var_time(row_nans,3));
                           PhTotalNT.Diastolic(i) = length(var_time(row_nans,3));
                        end
                    else % Other variables
%                         PhMsMeanPatient.(var_name)(i) = mean(cellfun(@(a) mean(a,'omitnan'),variables_list));
                        PhMsMedianPatient.(var_name)(i) = median(cellfun(@(a) median(a,'omitnan'),variables_list));
                        %PhMsMinPatient.(var_name)(i) = min(cell2mat(cellfun(@(a) min(a),variables_list,'UniformOutput',false)));
                        %PhMsMaxPatient.(var_name)(i) = max(cell2mat(cellfun(@(a) max(a),variables_list,'UniformOutput',false)));
                        PhMsLastPatient.(var_name)(i) = variables_list(end,1);
                        
                        % Check if measurement was taken during screening
                        in_screen = arrayfun(@(x) x >= max_screen_day & x <= min_screen_day , table2array(s(idx_conj,{'days_from_ok'})));
                        PhMsNMPatient.(var_name)(i) = sum(cellfun(@(a) ~any(isnan(a)),variables_list(in_screen)));
                        % Chech if measurement was taken before screening
                        % but up to a year
                        NofT = arrayfun(@(x) x >= -365 & x <= max_screen_day , table2array(s(idx_conj,{'days_from_ok'})));
                        PhNTPrior.(var_name)(i) = sum(cellfun(@(a) ~any(isnan(a)),variables_list(NofT)));
                        % Check if measurement was taken before a year of
                        % surgery
                        PhTotalNT.(var_name)(i) = sum(arrayfun(@(a) ~any(isnan(a)),table2array(s(idx_conj,{'days_from_ok'}))));
                    end
                else % No measurement for variable j
%                     PhMsMeanPatient.(var_name)(i) = NaN;
                    PhMsMedianPatient.(var_name)(i) = NaN;
                    %PhMsMinPatient.(var_name)(i) = NaN;
                    %PhMsMaxPatient.(var_name)(i) = NaN;
                    PhMsNMPatient.(var_name)(i) = 0;
                    PhMsLastPatient.(var_name)(i) = {NaN};
                    PhNTPrior.(var_name)(i) = 0;
                    PhTotalNT.(var_name)(i) = 0;
                end
            else % No measurement up to 1 year before surgery
%                 PhMsMeanPatient.(var_name)(i) = NaN;
                PhMsMedianPatient.(var_name)(i) = NaN;
                %PhMsMinPatient.(var_name)(i) = NaN;
                %PhMsMaxPatient.(var_name)(i) = NaN;
                PhMsNMPatient.(var_name)(i) = 0;
                PhMsLastPatient.(var_name)(i) = {NaN};
                PhNTPrior.(var_name)(i) = 0;
                PhTotalNT.(var_name)(i) = 0;
            end
        else % Patient doesn't has measurement
%             PhMsMeanPatient.(var_name)(i) = NaN;
            PhMsMedianPatient.(var_name)(i) = NaN;
            %PhMsMinPatient.(var_name)(i) = NaN;
            %PhMsMaxPatient.(var_name)(i) = NaN;
            PhMsNMPatient.(var_name)(i) = 0;
            PhMsLastPatient.(var_name)(i) = {NaN};
            PhNTPrior.(var_name)(i) = 0;
            PhTotalNT.(var_name)(i) = 0;
        end
    end
end

clearvars var_list prealloc_num vars timelapse var_time idx_var var_name ...
    s variables_list tbl2 tmp patient_ID n k j idx_patient_ID...
    idx_meas_below_a_times i id_patient_screen idx_conj NofT screen_day ...
    syst row_nans min_screen_day max_screen_day in_windowX in_twindow ...
    in_screen A dias

% Delete 'qID119058' Tension
% PhMsMeanPatient.qID119058 = [];
PhMsMedianPatient.qID119058 = [];
%PhMsMinPatient.qID119058 = [];
%PhMsMaxPatient.qID119058 = [];
PhMsNMPatient.qID119058 = [];
PhMsLastPatient.qID119058 = [];
PhNTPrior.qID119058 = [];
PhTotalNT.qID119058 = [];

%% Bar plot
% Total measurements done per variable

% Find how many questionIDs have each patient
f_zeros = sum(ismissing(PhTotalNT{:,2:end},0));
f_missing = sum(ismissing(PhTotalNT{:,2:end},NaN));
f_missing = sum(vertcat(f_missing,f_zeros));
tmp = PhTotalNT{:,2:end}; zero_idx = tmp <= 0;
tmp(zero_idx) = NaN;

f1=mean((PhTotalNT{:,2:end}),'omitnan');
f2=max(tmp);
f3=min(tmp);
f=vertcat(f3,f1,f2);

mean_bar = sscanf('E5CD00','%2x%2x%2x',[1 3])/255;
min_bar = sscanf('418A96','%2x%2x%2x',[1 3])/255;
max_bar = sscanf('D1002F','%2x%2x%2x',[1 3])/255;

f1_rgb = repmat(mean_bar,size(f1,2),1);
f2_rgb = repmat(max_bar,size(f2,2),1);
f3_rgb = repmat(min_bar,size(f3,2),1);

var_gradient = ["#E5CC00", "#E3BE03", "#E2B006", "#E1A309", "#DF950C",...
    "#DE880F", "#DD7A12", "#DB6C15", "#DA5F19", "#D9511C", "#D7441F",...
    "#D63622", "#D52825", "#D31B28", "#D20D2B", "#D1002F", "#CB0039",...
    "#C60043", "#C1004D", "#BC0058", "#B70062", "#B2006C", "#AD0077",...
    "#A70081", "#A2008B", "#9D0095", "#9800A0", "#9300AA", "#8E00B4",...
    "#8900BF"];
var_gradient_rgb = zeros(length(var_gradient),3);
    
for i = 1:length(var_gradient)
    color_str = var_gradient{i};
    color_bar = sscanf(color_str(2:end),'%2x%2x%2x',[1 3])/255;
    var_gradient_rgb(i,:) = color_bar;
end

figure('Name','Max Min Avg measurements per variable per patient');
g = bar(f','EdgeColor','none','FaceColor','flat');
g(1).CData = f3_rgb;
g(2).CData = f1_rgb;
g(3).CData = f2_rgb;
xlabel(['question ID or variable, n=',num2str(size(PhTotalNT{1,2:end},2))],...
    'FontSize',16,'FontWeight','bold');
ylabel('Observations','FontSize',16,'FontWeight','bold');
fig = gcf; ax = gca;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 28; %Width
fig.Position(4) = 20; %Height
ax.FontSize = 14;
legend({'Min/patient','Mean/patient','Max/patient'},'FontSize',16);
print([subfold,'/',d,'/Plots/PhMsNMS6'],'-deps')
print([subfold,'/',d,'/Plots/PhMsNMS6c'],'-depsc')
fig.PaperOrientation = 'landscape';
print([subfold,'/',d,'/Plots/PhMsNMS6'], '-dpdf', '-fillpage')


figure('Name','Missingness per variable');
b = bar(f_missing,'EdgeColor','none','FaceColor','flat','FaceAlpha',0.8);
b.CData = var_gradient_rgb;
hline = refline([0 200]);
hline.Color = 'w'; hline.LineWidth = 3;
xlabel(['question ID or variable, n=',num2str(size(PhTotalNT{1,2:end},2))],...
    'FontSize',16,'FontWeight','bold');
ylabel('Frequency','FontSize',16,'FontWeight','bold');
fig = gcf; ax = gca;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 28; %Width
fig.Position(4) = 20; %Height
ax.FontSize = 14;
print([subfold,'/',d,'/Plots/PhMsNMS6_nan'],'-deps')
print([subfold,'/',d,'/Plots/PhMsNMS6c_nan'],'-depsc')
fig.PaperOrientation = 'landscape';
print([subfold,'/',d,'/Plots/PhMsNMS6_nan'], '-dpdf', '-fillpage')

clear hline f_zeros tmp i color_bar color_str

%% Histogram of frequencies of measurement per variable

n_vars = size(PhMsNMPatient{1,2:end},2);

qIDS = PhMsNMPatient.Properties.VariableNames(2:end);

for p=1:n_vars
    figure('Name','Frequency of each Physical Measurement');
    z = histogram(PhTotalNT{:, p+1},...
        'BinMethod', 'integer', 'FaceAlpha',0.3);
    xlim([0 Inf])
    ylabel('Frequency')
    title(cellstr(qIDS(p)));
end

clearvars p;
%% Histogram per variable

% for p = 1:n_vars
%     figure('Name','Physical Measures (Median)');
%     try
%         q(p)= histogram(PhMsMedianPatient{:, p+1},...
%             'BinMethod', 'integer', 'FaceAlpha',0.3);
%     catch
%         q(p)= histogram(PhMsMedianPatient{:, p+1}, 'FaceAlpha',0.3);
%     end
%     % xlim([-1 Inf])
%     ylabel('Unit')
%     title(cellstr(qIDS(p)));
%     print([subfold,'/',d,'/Plots/PhMsS6Hist',char(cellstr(qIDS(p)))],'-deps')
%     print([subfold,'/',d,'/Plots/PhMsS6Hist',char(cellstr(qIDS(p))),'c'],'-depsc')
%     print([subfold,'/',d,'/Plots/PhMsS6Hist',char(cellstr(qIDS(p)))], '-dpdf', '-fillpage')
%     close all
% end

clearvars n_vars
%%
% _Created by Aldo Arévalo_ 