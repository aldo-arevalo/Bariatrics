%% Build Table with time-series variables
% Mean, median, last measurement before surgery, times measured, ...
% Maximum and minimum value

load([subfold,'/',d,'/ImportedTables.mat'],'total_patientsID', 'PPOSdata')

%% Delete BMI, HEIGHT, and WEIGHT variables
% These variables will be used to estimate pre-operative variables 

% Remove question IDs
    % Weight (gewicht)
    PPOSdata(ismissing(PPOSdata.questionID, 'vID148748'),:)=[];
    
    % BMI
    PPOSdata(ismissing(PPOSdata.questionID, 'vID148749'),:)=[];
    
    % Height (lengte)
    PPOSdata(ismissing(PPOSdata.questionID, 'vID148747'),:)=[];
    
%% Remove variables with less than 100 measurments/results

% Delete string/text type answers
%     PPOSdata.omschrijving01 = [];

% Extract questionID column
    %vrID = table2cell(PPOSdata(:,'questionID'));
    
    % Column^/variables names
%     qIDs = unique(PPOSdata.questionID);
%     tbl3={};
% 
%         for i=1:length(qIDs)
%             var = qIDs(i);
%             idx_vrID = strcmp(vrID, char(var));
% 
%             tbl3(i,1)= var;
%             tbl3(i,2)= num2cell(sum(idx_vrID));
%         end

    tbl3 = tabulate(PPOSdata.questionID);

    % Minimal frequency
    a=10;

    % Specify the height of the table and delete measurementes below 100 times
    idx_meas_below_a_times = [tbl3{:,2}]' < a;

    qIDs = tbl3(~idx_meas_below_a_times,1);
        
% Row's names
patients_list = total_patientsID;

clear i a var idx_vrID idx_meas_below_a_times total_patientsID

%% Pre-allocate matrix
prealloc_num = array2table(NaN(length(patients_list), length(qIDs)+1));

PPOSNMPatient = prealloc_num;
PPOSLastPatient = prealloc_num;
PPOSMedianPatient = prealloc_num;

% Add column names
var_list = vertcat('PatientCode',qIDs);

PPOSNMPatient.Properties.VariableNames = var_list;
PPOSLastPatient.Properties.VariableNames = var_list;
PPOSMedianPatient.Properties.VariableNames = var_list;

for i=1:length(patients_list)
    patient_ID = patients_list(i);
    idx_patient_ID = [PPOSdata.PatientCode] == patient_ID;
    
    PPOSNMPatient.PatientCode(i) = patient_ID;
    PPOSLastPatient.PatientCode(i) = patient_ID;
    PPOSMedianPatient.PatientCode(i) = patient_ID;
    
    s = PPOSdata(idx_patient_ID,:);
       
    for j = 1:length(qIDs)
        var_name = qIDs{j};
        idx_var = strcmp(s.questionID,var_name);
        
        if sum(idx_var)>= 1
                vars = table2array(s(idx_var,5:7));
                n = size(vars,1);
                variables_list = {};
                
                    for k = 1:n
                        tmp = vars(k,1:end-1);
                        tmp(isnan(tmp)) = [];
                        variables_list{k} = tmp;
                    end
                    
                    clearvars tmp n;
                    
                if isempty(variables_list{1,1}) == 0
                    
%                     PPOSMeanPatient.(var_name)(i) = mean(cellfun(@(a) mean(a,'omitnan'),variables_list));
%                     PPOSMinPatient.(var_name)(i) = min(cell2mat(cellfun(@(a) min(a),variables_list,'UniformOutput',false)));
%                     PPOSMaxPatient.(var_name)(i) = max(cell2mat(cellfun(@(a) max(a),variables_list,'UniformOutput',false)));
                    PPOSNMPatient.(var_name)(i) = sum(cellfun(@(a) ~any(isnan(a)),variables_list));
                    PPOSMedianPatient.(var_name)(i) = median(cellfun(@(a) median(a,'omitnan'),variables_list));
                    PPOSLastPatient.(var_name)(i) = cell2mat(variables_list(end,1));
   
                else
%                     PPOSMeanPatient.(var_name)(i) = NaN;
%                     PPOSMinPatient.(var_name)(i) = NaN;
%                     PPOSMaxPatient.(var_name)(i) = NaN;
                    PPOSNMPatient.(var_name)(i) = 0;
                    PPOSMedianPatient.(var_name)(i) = NaN;
                    PPOSLastPatient.(var_name)(i) = NaN;
                end
                    
        else
%             PPOSMeanPatient.(var_name)(i) = NaN;
%             PPOSMinPatient.(var_name)(i) = NaN;
%             PPOSMaxPatient.(var_name)(i) = NaN;
            PPOSNMPatient.(var_name)(i) = 0;
            PPOSMedianPatient.(var_name)(i) = NaN;
            PPOSLastPatient.(var_name)(i) = NaN;
        end
                    
    end
end

clearvars var_list prealloc_num vars timelapse var_time idx_var var_name ...
    j i k idx_patient_ID variables_list s vrID patient_ID

%% Bar plot
% Frequency of variable measurements

% Find how many questionIDs have each patient
f=sum(~ismissing(PPOSNMPatient{:,2:end},0));

var_gradient = ["#E5CC00", "#E1C407", "#DEBC0E", "#DAB415", "#D7AC1D",...
    "#D3A524", "#D09D2B", "#CC9533", "#C98D3A", "#C58541", "#C27E48", ...
    "#BE7650", "#BB6E57", "#B7665E", "#B45F66", "#B0586C", "#AD5172", ...
    "#AA4A79", "#A7437F", "#A43D85", "#A1368C", "#9E2F92", "#9B2898", ...
    "#98219F", "#951BA5", "#9214AB", "#8F0DB2", "#8C06B8", "#8900BF"];
var_gradient_rgb = zeros(length(var_gradient),3);
    
for i = 1:length(var_gradient)
    color_str = var_gradient{i};
    color_bar = sscanf(color_str(2:end),'%2x%2x%2x',[1 3])/255;
    var_gradient_rgb(i,:) = color_bar;
end

clearvars i color_bar color_str

figure('Name','Frequency of Measurement per patient');
g = bar(f,'EdgeColor','none','FaceColor','flat','FaceAlpha',0.8);
g.CData = var_gradient_rgb;
hline = refline([0 1000]);
hline.Color = 'w'; h.linewidth = 6;
xlabel(['Question ID or variables, n=',num2str(length(f))],...
    'FontSize',16,'FontWeight','bold');
ylabel('Frequency','FontSize',16,'FontWeight','bold');
fig = gcf; ax = gca;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 28; %Width
fig.Position(4) = 20; %Height
ax.FontSize = 14;
print([subfold,'/',d,'/Plots/FreqPPPOS'],'-deps')
print([subfold,'/',d,'/Plots/FreqPPOSc'],'-depsc')
fig.PaperOrientation = 'landscape';
print([subfold,'/',d,'/Plots/FrewPPOS'], '-dpdf','-bestfit')


clearvars hline

%% Histogram per variable
n_vars = length(f);

for p = 1:n_vars
    figure('Name','Pre-Surgical Screening (mean)');
    q(p)= histogram(PPOSLastPatient{:, p+1},...
        'BinMethod', 'integer', 'FaceAlpha',0.3);
    % xlim([-1 Inf])
    ylabel('Unit')
    title(cellstr(qIDs(p)));
end

%% Histogram of frequencies of measurement per variable

for p=1:n_vars
    figure('Name','Frequency of each Pre-Surgical Screening question');
    z = histogram(PPOSNMPatient{:, p+1},...
        'BinMethod', 'integer', 'FaceAlpha',0.3);
    xlim([-0.50 Inf])
    ylabel('Frequency')
    title(cellstr(qIDs(p)));
    
    print([subfold,'/',d,'/Plots/PPOS7Hist',char(cellstr(qIDs(p)))],'-deps')
    print([subfold,'/',d,'/Plots/PPOS7Hist',char(cellstr(qIDs(p))),'c'],'-depsc')
    print([subfold,'/',d,'/Plots/PPOS7Hist',char(cellstr(qIDs(p)))], '-dpdf', '-fillpage')
end

clearvars p n_vars f;

%%
% _Created by Aldo Arévalo_ 