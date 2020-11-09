%% Find Physical Measures per patient
% Objective: to make a unique list U with QOL codes which have a
% frequency >= a

A = exist('PhMsData','var');

if A == 0
    load([subfold,'/',d,'/ImportedTables.mat'],'PhMsData')
end

% Create valid column/variable names
PhMsData.STELLING = matlab.lang.makeValidName(PhMsData.STELLING);

% Obtain statistics from data gathered before surgery
tbl2=tabulate(PhMsData.STELLING);

% Specify the height of the table and delete measurementes below a times
idx_meas_below_a_times = [tbl2{:,2}]' < a;

U = tbl2(~idx_meas_below_a_times,1);

%% Initialize struct
% Extract Comments from Stelling
var_names = unique(U);

data_struct.meas = {};
data_struct.days_from_ok = NaN;

myStruct.PatientCode = NaN;
for i = 1:length(var_names)
    myStruct.(var_names{i}) = data_struct;
end

%% Fill in the structure for each patient

list_PatientID = unique(PhMsData.PatientCode);

PhMsPerPatient = repmat(myStruct,length(list_PatientID),1); % Pre-allocate structs for all patients

for i=1:length(list_PatientID)
    patient_ID = list_PatientID(i);
    
    PhMsPerPatient(i).PatientCode = patient_ID;
    
    idx_patient_ID = PhMsData.PatientCode == patient_ID;
    
    s = PhMsData(idx_patient_ID,:);
    
    var_names = unique(s.STELLING);
    
    for j = 1:length(var_names)
        var_name = var_names{j};
        idx_var = strcmp(s.STELLING,var_name);
        
        Ant = table2cell(s(idx_var,{'antwoord'}));
        Ant = cellfun(@str2double,Ant);
        varstmp = table2array(s(idx_var,{'value01','value02'}));
        vars = horzcat(varstmp,Ant);
        n = size(vars,1);
        
        variables_list = {};
        for k = 1:n
            variables_list{k} = vars(k,~isnan(vars(k,:)));
        end
        days_from_ok_list = s.days_from_ok(idx_var);
        
        PhMsPerPatient(i).(var_name).meas = variables_list;
        PhMsPerPatient(i).(var_name).days_from_ok = days_from_ok_list;
    end
     clear idx_var vars n s var_names idx_patient_ID var_name ...
         days_from_ok_list k j i data_struct variables_list patient_ID
end

field = tbl2(idx_meas_below_a_times,1);
PhMsPerPatient = rmfield(PhMsPerPatient,field);

clearvars idx_meas_below_a_times field A
%% Bar graph

C = sortrows(tbl2,2);
% Colormap
var_gradient = ["#E5CD00","#E2C804","#E0C309","#DEBE0D","#DCB912",...
    "#DAB416","#D8AF1B","#D5AA1F","#D3A524","#D1A128","#CF9C2D",...
    "#CD9732","#CB9236","#C88D3B","#C6883F","#C48344","#C27E48",...
    "#C07A4D","#BE7551","#BB7056","#B96B5A","#B7665F","#B56164",...
    "#B35C68","#B1576D","#AE5271","#AC4E76","#AA497A","#A8447F",...
    "#A63F83","#A43A88","#A1358C","#9F3091","#9D2B96","#9B279A",...
    "#99229F","#971DA3","#9418A8","#9213AC","#900EB1","#8E09B5",...
    "#8C04BA","#8A00BF"];

var_gradient_rgb = zeros(length(var_gradient),3);
    
for i = 1:length(var_gradient)
    color_str = var_gradient{i};
    color_bar = sscanf(color_str(2:end),'%2x%2x%2x',[1 3])/255;
    var_gradient_rgb(i,:) = color_bar;
end

clearvars i color_bar color_str

figure('Name','Frequency of Physical Measures')
b=bar(cell2mat(C(:,2)),'EdgeColor','none','FaceColor','flat');
b.CData = var_gradient_rgb;
xlabel('Physical Measure','FontSize',16,'FontWeight','bold');
ylabel('Frequency','FontSize',16,'FontWeight','bold');
fig = gcf; ax = gca;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 28; %Width
fig.Position(4) = 20; %Height
ax.FontSize = 14;
print([subfold,'/',d,'/Plots/PhMeasBarsS3'],'-deps')
print([subfold,'/',d,'/Plots/PhMeasBarsS3c'],'-depsc')
fig.PaperOrientation = 'landscape';
print([subfold,'/',d,'/Plots/PhMeasBarsS3'], '-dpdf', '-fillpage')

clearvars C U x;