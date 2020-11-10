%% Build Table with mean values

A = exist('total_patientsID','var');

if A == 0
    load([subfold,'/',d,'/ImportedTables.mat'], 'total_patientsID')
end

A = exist('LabResultsPerPatient','var');

if A == 0
    load([subfold,'/',d,'/LabResultsStructure.mat'], 'LabResultsPerPatient')
end

clearvars A ASAT BKR00b Albumin ALAT ReferenceValues Calcium i
%% Obtain screening date per patient
lab_pack = {'PatientCode','BHE001','BHE000','BER002','BMC000','BMC002',...
    'BTR006','BLE001','BGL003','BBI001','BAS002','BAL015','BLD004','BAL014',...
    'BGT001','BUR002','BKR000b','BKR015','BKA000','BNA001','BCA006b',...
    'BFO002','BAL002','BCR002','BCH004','BHD000','BCH014','BLD002',...
    'BPR012','BHE030','BIN000','BCP002','BPA016','BTS003','BCO017',...
    'BIJ001','BFE006','BFO000','BMA002','ScreenDay'};

% Rows
u = total_patientsID; % All patients

LabPackTimesPatient = array2table(NaN(length(u), length(lab_pack)));

LabPackTimesPatient.Properties.VariableNames = lab_pack;

for i=1:length(u)
    patient_ID = u(i);
    idx_patient_ID = [LabResultsPerPatient.PatientCode] == patient_ID;
    
    LabPackTimesPatient.PatientCode(i) = patient_ID;
    
    for j=1:(length(lab_pack)-2)
        var_name = lab_pack{j+1};
        
        A = [LabResultsPerPatient(idx_patient_ID).(var_name)];
        
        if isempty(A) == 0
           n_values = size(A.Result,1);
           unique_days = unique(A.days_from_ok);
           
           switch (size(unique_days,1))
               case 1
                   vidx = arrayfun(@(x) x >= -365, unique_days);
                   if sum(vidx) >= 1
                       LabPackTimesPatient.(var_name)(i) = unique_days(vidx,1);
                   else
                       LabPackTimesPatient.(var_name)(i) = NaN;
                   end
               case 2
                   vidx = arrayfun(@(x) x >= -365, unique_days);
                   if sum(vidx) >= 1
                       mat_tmp = unique_days(vidx);
                       LabPackTimesPatient.(var_name)(i) = mat_tmp(end);
                       clearvars mat_tmp
                   else
                      LabPackTimesPatient.(var_name)(i) = NaN;
                   end
               otherwise % >= 3
                   vidx = arrayfun(@(x) x >= -365, unique_days);
                   if sum(vidx) >= 1
                       mat_tmp = unique_days(vidx);
                       LabPackTimesPatient.(var_name)(i) = mat_tmp(end);
                       clearvars mat_tmp
                   else
                      LabPackTimesPatient.(var_name)(i) = NaN;
                   end
           end
            clearvars vidx
        else % empty == 1
            LabPackTimesPatient.(var_name)(i) = NaN;
        end
    end
    
    % Find the screening day
    LabPackTimesPatient.ScreenDay(i) = median(table2array(LabPackTimesPatient(i,2:end-1)),'omitnan');
end

clearvars i j patientidx patientID n_values unique_days var_name
%% Columns
var_list = fieldnames(LabResultsPerPatient);

    % Remove variables of sCr and total serum calcium
    % var_list(strcmp(var_list, 'BKR000'),:)=[]; % Creatinine
    var_list(strcmp(var_list, 'BCA006'),:)=[]; % Calcium (quelated)
    var_list(strcmp(var_list, 'BAM012'),:)=[]; % amylase
    var_list(strcmp(var_list, 'BAN037'),:)=[]; % 
    var_list(strcmp(var_list, 'BAP012'),:)=[]; % APTT
    var_list(strcmp(var_list, 'BBA002'),:)=[]; % basofielen
    var_list(strcmp(var_list, 'BBA004'),:)=[]; % base excess
    var_list(strcmp(var_list, 'BBE001'),:)=[]; % bezinking erytrocyten
    var_list(strcmp(var_list, 'BBI000'),:)=[]; % bilirubine geconjugeerd
    var_list(strcmp(var_list, 'BBI002'),:)=[]; % bicarbonaat
    var_list(strcmp(var_list, 'BCH003'),:)=[]; % chloride
    var_list(strcmp(var_list, 'BCK003'),:)=[]; % CK-MB
    var_list(strcmp(var_list, 'BCK004'),:)=[]; % CK
    var_list(strcmp(var_list, 'BCO017'),:)=[]; % plasma cortisol
    var_list(strcmp(var_list, 'BCO029'),:)=[]; % totaal CO2
    var_list(strcmp(var_list, 'BDI003'),:)=[]; % differentiele telling
    var_list(strcmp(var_list, 'BEI000'),:)=[]; % eiwit totaal
    var_list(strcmp(var_list, 'BEO001'),:)=[]; % eosinofielen
    var_list(strcmp(var_list, 'BFI005'),:)=[]; % fibrinogeen
    var_list(strcmp(var_list, 'BFI007'),:)=[]; % FIO2
    var_list(strcmp(var_list, 'BFS002'),:)=[]; % FSH
    var_list(strcmp(var_list, 'BHC004'),:)=[]; % hCG (+ beta-subunits)
    var_list(strcmp(var_list, 'BHE016'),:)=[]; % hemoglobine A1c (NGSP)
    var_list(strcmp(var_list, 'BHO001'),:)=[]; % homocysteine
    var_list(strcmp(var_list, 'BIN001'),:)=[]; % inhalatie-allergie screening
    var_list(strcmp(var_list, 'BIN008'),:)=[]; % INR (PT)
    var_list(strcmp(var_list, 'BKR011'),:)=[]; % kruisbloed
    var_list(strcmp(var_list, 'BKR015'),:)=[]; % kreatinine klaring MDRD
    var_list(strcmp(var_list, 'BLH002'),:)=[]; % LH
    var_list(strcmp(var_list, 'BLOEDGRR'),:)=[]; % bloedgroep rhesus D conclusie
    var_list(strcmp(var_list, 'BLY003'),:)=[]; % lymfocyten
    var_list(strcmp(var_list, 'BMO004'),:)=[]; % monocyten
    var_list(strcmp(var_list, 'BNE000'),:)=[]; % neutrofielen
    var_list(strcmp(var_list, 'BOE001'),:)=[]; % oestradiol
    var_list(strcmp(var_list, 'BPC000'),:)=[]; % pCO2
    var_list(strcmp(var_list, 'BPH000'),:)=[]; % pH
    var_list(strcmp(var_list, 'BPO000'),:)=[]; % pO2
    var_list(strcmp(var_list, 'BPR012'),:)=[]; % protrombine tijd
    var_list(strcmp(var_list, 'BRE004'),:)=[]; % reticulocyten
    var_list(strcmp(var_list, 'BRH003'),:)=[]; % rhesus D
    var_list(strcmp(var_list, 'BSE001'),:)=[]; % serologie (PAMM)
    var_list(strcmp(var_list, 'BSO001'),:)=[]; % sO2
    var_list(strcmp(var_list, 'BST011'),:)=[]; % standaard bicarbonaat
    var_list(strcmp(var_list, 'BTE002'),:)=[]; % temperatuur
    var_list(strcmp(var_list, 'BTR001'),:)=[]; % transferrine
    var_list(strcmp(var_list, 'BTR021'),:)=[]; % transferrineverzadiging
    var_list(strcmp(var_list, 'BTR031'),:)=[]; % c-troponine-I
    var_list(strcmp(var_list, 'BUR001'),:)=[]; % uraat
    var_list(strcmp(var_list, 'OBA000'),:)=[]; % 
    var_list(strcmp(var_list, 'OBI000'),:)=[]; % 
    var_list(strcmp(var_list, 'OBL000'),:)=[]; %
    var_list(strcmp(var_list, 'OCA000'),:)=[]; %
    var_list(strcmp(var_list, 'OCA001'),:)=[]; %
    var_list(strcmp(var_list, 'OCH000'),:)=[]; %
    var_list(strcmp(var_list, 'OGL000'),:)=[]; %
    var_list(strcmp(var_list, 'OHE000'),:)=[]; %
    var_list(strcmp(var_list, 'OKA000'),:)=[]; %
    var_list(strcmp(var_list, 'OLA000'),:)=[]; %
    var_list(strcmp(var_list, 'OLI001'),:)=[]; %
    var_list(strcmp(var_list, 'OME000'),:)=[]; %
    var_list(strcmp(var_list, 'ONA000'),:)=[]; %
    var_list(strcmp(var_list, 'OPC000'),:)=[]; %
    var_list(strcmp(var_list, 'OPH000'),:)=[]; %
    var_list(strcmp(var_list, 'OPO000'),:)=[]; %
    var_list(strcmp(var_list, 'OSO000'),:)=[]; %
    var_list(strcmp(var_list, 'OST000'),:)=[]; %
    var_list(strcmp(var_list, 'OTE000'),:)=[]; %
    var_list(strcmp(var_list, 'OVO000'),:)=[]; %
    var_list(strcmp(var_list, 'TGL000'),:)=[]; %
    var_list(strcmp(var_list, 'UAL001'),:)=[]; %
    var_list(strcmp(var_list, 'UAL002'),:)=[]; %
    var_list(strcmp(var_list, 'UAL004'),:)=[]; %
    var_list(strcmp(var_list, 'UAL005'),:)=[]; %
    var_list(strcmp(var_list, 'UBA000'),:)=[]; %
    var_list(strcmp(var_list, 'UBI000'),:)=[]; %
    var_list(strcmp(var_list, 'UEI001'),:)=[]; %
    var_list(strcmp(var_list, 'UEI002'),:)=[]; %
    var_list(strcmp(var_list, 'UEI003'),:)=[]; %
    var_list(strcmp(var_list, 'UEI006'),:)=[]; %
    var_list(strcmp(var_list, 'UER000'),:)=[]; %
    var_list(strcmp(var_list, 'UER002'),:)=[]; %
    var_list(strcmp(var_list, 'UGL000'),:)=[]; %
    var_list(strcmp(var_list, 'UKE000'),:)=[]; %
    var_list(strcmp(var_list, 'UKR002'),:)=[]; %
    var_list(strcmp(var_list, 'UKR004'),:)=[]; %
    var_list(strcmp(var_list, 'ULE000'),:)=[]; %
    var_list(strcmp(var_list, 'UNA002'),:)=[]; %
    var_list(strcmp(var_list, 'UNI000'),:)=[]; %
    var_list(strcmp(var_list, 'UPH000'),:)=[]; %
    var_list(strcmp(var_list, 'UPL000'),:)=[]; %
    var_list(strcmp(var_list, 'USE000'),:)=[]; %
    var_list(strcmp(var_list, 'USO000'),:)=[]; %
    var_list(strcmp(var_list, 'USO001'),:)=[]; %
    var_list(strcmp(var_list, 'UUR007'),:)=[]; %
    var_list(strcmp(var_list, 'UVO000'),:)=[]; %
    var_list(strcmp(var_list, 'ASAT_ALAT'),:)=[]; %

    % Lab package list
    lab_pack = {'BHE001','BHE000','BER002','BMC000','BMC002',...
    'BTR006','BLE001','BGL003','BBI001','BAS002','BAL015','BLD004','BAL014',...
    'BGT001','BUR002','BKR000b','BKA000','BNA001','BCA006b',...
    'BFO002','BAL002','BCR002','BCH004','BHD000','BCH014','BLD002',...
    'BHE030','BIN000','BCP002','BPA016','BTS003',...
    'BIJ001','BFE006','BFO000','BMA002'};
    
% Pre-allocate matrix
prealloc_num = array2table(NaN(length(u), length(var_list)));

LabMeanPatient = prealloc_num;
LabMinPatient = prealloc_num;
LabMaxPatient = prealloc_num;
LabNMPatient = prealloc_num;
LabLastPatient = prealloc_num;
LabMedianPatient = prealloc_num;
LabNTPrior = prealloc_num;
LabTotalNT = prealloc_num;

ScreenTimeMax = prealloc_num;
ScreenTimeMedian = prealloc_num;
ScreenTimeLast = prealloc_num;

% Add column names
LabMeanPatient.Properties.VariableNames = var_list;
LabMinPatient.Properties.VariableNames = var_list;
LabMaxPatient.Properties.VariableNames = var_list;
LabNMPatient.Properties.VariableNames = var_list;
LabLastPatient.Properties.VariableNames = var_list;
LabMedianPatient.Properties.VariableNames = var_list;
LabNTPrior.Properties.VariableNames = var_list;
LabTotalNT.Properties.VariableNames = var_list;

ScreenTimeMax.Properties.VariableNames = var_list;
ScreenTimeMedian.Properties.VariableNames = var_list;
ScreenTimeLast.Properties.VariableNames = var_list;

No_features = length(var_list);

for i=1:length(u)
    patient_ID = u(i);
    idx_patient_ID = [LabResultsPerPatient.PatientCode] == patient_ID;
    id_patient_screen = table2array(LabPackTimesPatient(:,{'PatientCode'})) == patient_ID;
    screen_day = table2array(LabPackTimesPatient(id_patient_screen,{'ScreenDay'}));
    max_screen_day = screen_day - 15;
    min_screen_day = screen_day + 15;
    
    LabMeanPatient.PatientCode(i) = patient_ID;
    LabMinPatient.PatientCode(i) = patient_ID;
    LabMaxPatient.PatientCode(i) = patient_ID;
    LabNMPatient.PatientCode(i) = patient_ID;
    LabLastPatient.PatientCode(i) = patient_ID;
    LabMedianPatient.PatientCode(i) = patient_ID;
    LabNTPrior.PatientCode(i) = patient_ID;
    LabTotalNT.PatientCode(i) = patient_ID;
    
    ScreenTimeMax.PatientCode(i) = patient_ID;
    ScreenTimeMedian.PatientCode(i) = patient_ID;
    ScreenTimeLast.PatientCode(i) = patient_ID;
    
    for j = 1:(length(var_list)-1)
        var_name = var_list{j+1};
        A = [LabResultsPerPatient(idx_patient_ID).(var_name)];
        
        if isempty(A) == 0
            if isempty(A.Result) == 0
                % Find indexes of values that are in the timewindow of
                % screening               
                in_twindow = arrayfun(@(x) x >= -365, A.days_from_ok);
                
                if sum(in_twindow) >= 1
                    in_screen = arrayfun(@(x) x >= max_screen_day & x <= min_screen_day , A.days_from_ok);
                    if sum(in_screen) >= 1
                        ScreenTimeMax.(var_name)(i) = max(A.days_from_ok(in_screen,1)); %index of closest value to zero
                        ScreenTimeLast.(var_name)(i) = min(A.days_from_ok(in_screen,1));
                        ScreenTimeMedian.(var_name)(i) = median(A.days_from_ok(in_screen,1),'omitnan');

                        LabMeanPatient.(var_name)(i) = mean(A.Result(in_screen,1),'omitnan');
                        LabMinPatient.(var_name)(i) = min(A.Result(in_screen,1));
                        LabMaxPatient.(var_name)(i) = max(A.Result(in_screen,1));
                        LabMedianPatient.(var_name)(i) = median(A.Result(in_screen,1),'omitnan');
                        LabLastPatient.(var_name)(i) = A.Result(find(in_screen, 1, 'last' ),1);
                        LabNMPatient.(var_name)(i) = sum(~isnan(A.Result(in_screen,1)));
                    else
                        LabMeanPatient.(var_name)(i) = NaN;
                        LabMinPatient.(var_name)(i) = NaN;
                        LabMaxPatient.(var_name)(i) = NaN;
                        LabMedianPatient.(var_name)(i) = NaN;
                        LabLastPatient.(var_name)(i) = NaN;
                        ScreenTimeMax.(var_name)(i) = NaN;
                        ScreenTimeMedian.(var_name)(i) = NaN;
                        ScreenTimeLast.(var_name)(i) = NaN;
                    end
                    
                    % Find values measured before screening data and up to
                    % a year before surgery
                    NofT = arrayfun(@(x) x >= -365 & x <= max_screen_day , A.days_from_ok);
                    LabNTPrior.(var_name)(i) = sum(~isnan(A.Result(NofT,1)));
                    LabTotalNT.(var_name)(i) = sum(~isnan(A.Result(in_twindow,1)));
                else
                    LabMeanPatient.(var_name)(i) = NaN;
                    LabMinPatient.(var_name)(i) = NaN;
                    LabMaxPatient.(var_name)(i) = NaN;
                    LabNMPatient.(var_name)(i) = 0;
                    LabMedianPatient.(var_name)(i) = NaN;
                    LabLastPatient.(var_name)(i) = NaN;
                    LabNTPrior.(var_name)(i) = 0;
                    LabTotalNT.(var_name)(i) = 0;

                    ScreenTimeMax.(var_name)(i) = NaN;
                    ScreenTimeMedian.(var_name)(i) = NaN;
                    ScreenTimeLast.(var_name)(i) = NaN;
                end
            else
                LabMeanPatient.(var_name)(i) = NaN;
                LabMinPatient.(var_name)(i) = NaN;
                LabMaxPatient.(var_name)(i) = NaN;
                LabNMPatient.(var_name)(i) = 0;
                LabMedianPatient.(var_name)(i) = NaN;
                LabLastPatient.(var_name)(i) = NaN;
                LabNTPrior.(var_name)(i) = 0;
                LabTotalNT.(var_name)(i) = 0;
                
                ScreenTimeMax.(var_name)(i) = NaN;
                ScreenTimeMedian.(var_name)(i) = NaN;
                ScreenTimeLast.(var_name)(i) = NaN;
            end
        else
            LabMeanPatient.(var_name)(i) = NaN;
            LabMinPatient.(var_name)(i) = NaN;
            LabMaxPatient.(var_name)(i) = NaN;
            LabNMPatient.(var_name)(i) = 0;
            LabMedianPatient.(var_name)(i) = NaN;
            LabLastPatient.(var_name)(i) = NaN;
            LabNTPrior.(var_name)(i) = 0;
            LabTotalNT.(var_name)(i) = 0;
            
            ScreenTimeMax.(var_name)(i) = NaN;
            ScreenTimeMedian.(var_name)(i) = NaN;
            ScreenTimeLast.(var_name)(i) = NaN;
        end
    end
   clearvars var_name idx_closeTOzero A idx_patient_ID in_twindow NofT ...
       screen_day;
end

clearvars prealloc_num j i patient_ID id_patient_screen ...
    max_screen_day min_screen_day in_screen;

%% Number of additional measurements for each parameter per patient

AddLabScreen = LabNMPatient;
AddLabScreen(:,lab_pack) =[];

for i=1:size(AddLabScreen,2)-1
    idx_greater2 = table2array(AddLabScreen(:,i+1))>=2;
    if sum(idx_greater2) >=1
        AddLabScreen(idx_greater2,i+1) = array2table(1);
        AddLabScreen(~idx_greater2,i+1) = array2table(0);
    else
        AddLabScreen(:,i+1) = array2table(0);
    end
end

AddScore = array2table(sum(table2array(AddLabScreen(:,2:end)),2));
AddLabScreen = horzcat(AddLabScreen,AddScore);
AddLabScreen.Properties.VariableNames{end} = 'AddScore';
AddScore = array2table(sum(table2array(LabNMPatient(:,2:end)),2));
AddLabScreen = horzcat(AddLabScreen,AddScore);
AddLabScreen.Properties.VariableNames{end} = 'SumOFAdd';

AddLabBeforeScreen = LabNTPrior;
AddLabBeforeScreen(:,lab_pack) =[];

for i=1:size(AddLabBeforeScreen,2)-1
    idx_greater2 = table2array(AddLabBeforeScreen(:,i+1))>=2;
    if sum(idx_greater2) >=1
        AddLabBeforeScreen(idx_greater2,i+1) = array2table(1);
        AddLabBeforeScreen(~idx_greater2,i+1) = array2table(0);
    else
        AddLabBeforeScreen(:,i+1) = array2table(0);
    end
end

AddScore = array2table(sum(table2array(AddLabBeforeScreen(:,2:end)),2));
AddLabBeforeScreen = horzcat(AddLabBeforeScreen,AddScore);
AddLabBeforeScreen.Properties.VariableNames{end} = 'AddScore';
AddScore = array2table(sum(table2array(LabNTPrior(:,2:end)),2));
AddLabBeforeScreen = horzcat(AddLabBeforeScreen,AddScore);
AddLabBeforeScreen.Properties.VariableNames{end} = 'SumOFAdd';

AddLab1year = LabTotalNT;
AddLab1year(:,lab_pack) =[];

for i=1:size(AddLab1year,2)-1
    idx_greater2 = table2array(AddLab1year(:,i+1))>=2;
    if sum(idx_greater2) >=1
        AddLab1year(idx_greater2,i+1) = array2table(1);
        AddLab1year(~idx_greater2,i+1) = array2table(0);
    else
        AddLab1year(:,i+1) = array2table(0);
    end
end

AddScore = array2table(sum(table2array(AddLab1year(:,2:end)),2));
AddLab1year = horzcat(AddLab1year,AddScore);
AddLab1year.Properties.VariableNames{end} = 'AddScore';
AddScore = array2table(sum(table2array(LabTotalNT(:,2:end)),2));
AddLab1year = horzcat(AddLab1year,AddScore);
AddLab1year.Properties.VariableNames{end} = 'SumOFAdd';

clearvars i idx_greater2 AddScore
%% Bar plot
% Frequency of variable measurements

% Find how many Lab analyses have each patient
f=sum(~ismissing(LabNMPatient{:,2:end},0)); % Lab frequency

% 131 colors
% var_gradient = ["#E5CC00", "#E4C800", "#E4C501", "#E4C201", "#E3BF02",...
%     "#E3BC03", "#E3B903", "#E2B604", "#E2B204", "#E2AF05", "#E2AC06",...
%     "#E1A906", "#E1A607", "#E1A308", "#E0A008", "#E09C09", "#E09909",...
%     "#E0960A", "#DF930B", "#DF900B", "#DF8D0C", "#DE8A0C", "#DE860D",...
%     "#DE830E", "#DD800E", "#DD7D0F", "#DD7A10", "#DD7710", "#DC7411",...
%     "#DC7011", "#DC6D12", "#DB6A13", "#DB6713", "#DB6414", "#DB6114",...
%     "#DA5E15", "#DA5B16", "#DA5716", "#D95417", "#D95118", "#D94E18",...
%     "#D94B19", "#D84819", "#D8451A", "#D8411B", "#D73E1B", "#D73B1C",...
%     "#D7381C", "#D6351D", "#D6321E", "#D62F1E", "#D62B1F", "#D52820",...
%     "#D52520", "#D52221", "#D41F21", "#D41C22", "#D41923", "#D41523",...
%     "#D31224", "#D30F24", "#D30C25", "#D20926", "#D20626", "#D20327",...
%     "#D20028", "#D0002A", "#CF002C", "#CE002E", "#CD0031", "#CC0033",...
%     "#CB0035", "#CA0038", "#C9003A", "#C7003C", "#C6003F", "#C50041",...
%     "#C40043", "#C30046", "#C20048", "#C1004A", "#C0004D", "#BE004F",...
%     "#BD0051", "#BC0054", "#BB0056", "#BA0058", "#B9005B", "#B8005D",...
%     "#B7005F", "#B50062", "#B40064", "#B30066", "#B20069", "#B1006B",...
%     "#B0006D", "#AF0070", "#AE0072", "#AC0074", "#AB0076", "#AA0079",...
%     "#A9007B", "#A8007D", "#A70080", "#A60082", "#A50084", "#A30087",...
%     "#A20089", "#A1008B", "#A0008E", "#9F0090", "#9E0092", "#9D0095",...
%     "#9C0097", "#9A0099", "#99009C", "#98009E", "#9700A0", "#9600A3",...
%     "#9500A5", "#9400A7", "#9300AA", "#9100AC", "#9000AE", "#8F00B1",...
%     "#8E00B3", "#8D00B5", "#8C00B8", "#8B00BA", "#8A00BC", "#8900BF"];

% 44 colors

var_gradient = {'#cb356b','#ca3569','#ca3568','#ca3567','#c93565',...
    '#c93664','#c93663','#c83661','#c83660','#c8375f','#c7375d',...
    '#c7375c','#c7375b','#c63859','#c63858','#c63857','#c53855',...
    '#c53854','#c53953','#c43951','#c43950','#c4394f','#c33a4d',...
    '#c33a4c','#c33a4b','#c23a49','#c23b48','#c23b47','#c13b45',...
    '#c13b44','#c13b43','#c03c41','#c03c40','#c03c3f','#bf3c3d',...
    '#bf3d3c','#bf3d3b','#be3d39','#be3d38','#be3e37','#bd3e35',...
    '#bd3e34','#bd3e33'};

var_gradient_rgb = zeros(length(var_gradient),3);
    
for i = 1:length(var_gradient)
    color_str = var_gradient{i};
    color_bar = sscanf(color_str(2:end),'%2x%2x%2x',[1 3])/255;
    var_gradient_rgb(i,:) = color_bar;
end

clearvars i color_bar color_str

figure('Name','Frequency of Lab Analyses per patient');
g = bar(sort(f),'EdgeColor','none','FaceColor','flat');
g.CData = var_gradient_rgb;
hline = refline([0 1000]);
hline.Color = 'w'; hline.LineWidth = 3;
xlabel(['Lab analyses, n=',num2str(length(f))],'FontSize',16,...
    'FontWeight','bold');
ylabel('Frequency','FontSize',16,'FontWeight','bold');
fig = gcf; ax = gca;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 28; %Width
fig.Position(4) = 20; %Height
ax.FontSize = 14;
print([subfold,'/',d,'/Plots/LabS5'],'-deps')
print([subfold,'/',d,'/Plots/LabS5c'],'-depsc')
fig.PaperOrientation = 'landscape';
print([subfold,'/',d,'/Plots/LabS5'], '-dpdf', '-fillpage')

%% Histogram per variable
n_vars = length(f);
var_list = LabMeanPatient.Properties.VariableNames;

for p = 1:n_vars
%     figure('Name','Lab Results (median)');
%     q2(p)= histogram(LabMedianPatient{:, p+1},...
%         'BinMethod', 'integer', 'FaceAlpha',0.3);
%     ylabel('Unit')
%     title(cellstr(var_list(p+1)));
%     print([subfold,'/',d,'/Plots/',char(cellstr(var_list(p+1))),'median'],'-deps')
%     print([subfold,'/',d,'/Plots/',char(cellstr(var_list(p+1))),'medianc'],'-depsc')
%     print([subfold,'/',d,'/Plots/',char(cellstr(var_list(p+1))),'median'], '-dpdf', '-fillpage')
    
    figure('Name','Lab Results (Last measurements)');
    q3(p)= histogram(LabLastPatient{:, p+1},...
        'BinMethod', 'integer', 'FaceAlpha',0.3);
    ylabel('Unit')
    title(cellstr(var_list(p+1)));
    print([subfold,'/',d,'/Plots/',char(cellstr(var_list(p+1))),'_last'],'-deps')
    print([subfold,'/',d,'/Plots/',char(cellstr(var_list(p+1))),'_lastc'],'-depsc')
    print([subfold,'/',d,'/Plots/',char(cellstr(var_list(p+1))),'_last'], '-dpdf', '-fillpage')

    close all
end

clearvars n_vars f p q q2 q3
%% Lab package variables

for p=1:(length(lab_pack)-2)
    var_name = lab_pack{p+1};
    
    figure('Name',['Screening times ',var_name]);
    vidx = arrayfun(@(x) x >= -365 & x < 1, ScreenTimeLast{:, var_name});
    c=-365:12:0;      
    q4(p)= histogram(ScreenTimeLast{vidx, var_name},c,'DisplayStyle','stairs');
    hold on
    
    vidx = arrayfun(@(x) x >= -365 & x < 1, ScreenTimeMax{:, var_name});
    q6(p)= histogram(ScreenTimeMax{vidx, var_name},c,'DisplayStyle','stairs');
    
    vidx = arrayfun(@(x) x >= -365 & x < 1, ScreenTimeMedian{:, var_name});
    q8(p)= histogram(ScreenTimeMedian{vidx, var_name},c,'DisplayStyle','stairs');
    hold off
    
    ylabel('Patients')
    xlabel('Time (days)')
    title(cellstr(var_list(p+1)));
    
    legend('Last','Max','Median','Location','Best')
    
    print([subfold,'/',d,'/Plots/',char(cellstr(var_list(p+1))),'screen'],'-deps')
    print([subfold,'/',d,'/Plots/',char(cellstr(var_list(p+1))),'screenc'],'-depsc')
    print([subfold,'/',d,'/Plots/',char(cellstr(var_list(p+1))),'screen'], '-dpdf', '-fillpage')
    
    close all
end

clearvars q4 q6 q8 p c hline g PhMsPerPatient var_name vidx i

%%
% _% Created by Aldo Arévalo_ 