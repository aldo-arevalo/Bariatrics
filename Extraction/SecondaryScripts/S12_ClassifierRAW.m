%% Code to estimate Sucess or no-sucess class
% Sources: Pre-operativeValues2018 and FollowUpDataRAW

load([subfold,'/',d,'/PreOperativeValuesRAW.mat'],'weight_max')
load([subfold,'/',d,'/RawData/FollowUpDataRAW.mat'], 'FollowUpData')
load([subfold,'/',d,'/RawData/GeneralDataRAW.mat'])

% Convert to numeric column days_from_ok
Ant = table2cell(FollowUpData(:,{'days_from_ok'}));
FollowUpData(:,{'days_from_ok'}) = [];
FollowUpData(:,{'days_from_ok'}) = array2table(cellfun(@str2double,Ant));

clearvars Ant

% Row's names
patients_list = unique(GeneralData.PatientNr);

% Column's names
var_names = {'TWL','WL','kg'};

%% Pre-allocate matrix
preallocate_nan = array2table(NaN(length(patients_list), ...
    length(var_names)+1));

Classifiers1yr = preallocate_nan;
Classifiers2yr = preallocate_nan;

% Add column names
col_names = horzcat('PatientCode',var_names);
Classifiers1yr.Properties.VariableNames = col_names;
Classifiers2yr.Properties.VariableNames = col_names;

clearvars preallocate_nan col_names

for i=1:length(patients_list)
    patient_ID = patients_list(i);
    
    Classifiers1yr.PatientCode(i) = patient_ID;
    Classifiers2yr.PatientCode(i) = patient_ID;
    
    
    % Pre-Operative values
    idx_patientPRE = [weight_max.PatientCode] == patient_ID;
    Weight_patientPRE = weight_max.Max(idx_patientPRE);
    
    
    % Follow-Up values
    idx_patientFP = [FollowUpData.PatientCode] == patient_ID;
    
        % Weight
        Ws_FP = FollowUpData(idx_patientFP,{'folgewicht','days_from_ok'});
        
        if isempty(Ws_FP) == 0
            n = size(Ws_FP,1);
                weights_FP = {};
                for k = 1:n
                    tmp = table2array(Ws_FP(k,1:end-1));
                    tmp(isnan(tmp)) = [];
                    weights_FP{k} = tmp;
                end

            weights_patients = horzcat( weights_FP',table2cell(Ws_FP(:,end)));
            weights_patients(any(cellfun(@isempty,weights_patients),2),:) = [];
            weights_patients(cellfun(@(x) x < 41.3817 | x > 300, weights_patients(:,1)),:) = [];
            weights_patients = sortrows(weights_patients,2);
            clearvars tmp n k;

            % 1st year >=335 days & <395 days
            % 1st year >=180 days (6 months) & <545 days (1 yr + 6 months)
            % 1st year >=270 days (9 months) & <450 days (1 yr + 3 months)
            % Minimum value
            vidx = cellfun(@(x) x >= 270 & x < 450, weights_patients(:,2));
                if sum(vidx) == 0
                    weight_1yr = NaN;
                else
                    weight_1yr = min(cell2mat(weights_patients(vidx,1)));
                end

            % 2nd year >= 700 days & < 760 days
            % 2nd year >= 545 days & < 910 days (2 yr + 6 months)
            % 2nd year >= 630 days & < 810 days (2 yr + 3 months)
            % Minimum value
            vidx2 = cellfun(@(x) x >= 630 & x <810, weights_patients(:,2));

                if sum(vidx2) == 0
                    weight_2yr = NaN;
                else
                    weight_2yr = min(cell2mat(weights_patients(vidx2,1)));
                end
        else
            weight_1yr = NaN;
            weight_2yr = NaN;
        end
            clearvars vidx2 vidx weights_FP Ws_FP
        
            
        
            
   Classifiers1yr.TWL(i) = (Weight_patientPRE - weight_1yr)/Weight_patientPRE;
   Classifiers1yr.WL(i) = Weight_patientPRE - weight_1yr;
   Classifiers1yr.kg(i) = weight_1yr;
   
   Classifiers2yr.TWL(i) = (Weight_patientPRE - weight_2yr)/Weight_patientPRE;
   Classifiers2yr.WL(i) = Weight_patientPRE - weight_2yr;
   Classifiers2yr.kg(i) = weight_2yr;
   
clearvars Weight_patientPRE weight_1yr weight_2yr patient_ID idx_patientFP ...
    idx_patientPRE weights_patients
end

% Delete negative values
% Classifiers2yr((Classifiers2yr.TWL<=0),:)=[];
% Classifiers1yr((Classifiers1yr.TWL<=0),:)=[];

% Convert negative values to zero
varfun(@(x) x < 0, Classifiers1yr(:,{'TWL'}))
varfun(@(x) x < 0, Classifiers2yr(:,{'TWL'}))

save([subfold,'/',d,'/Classifiers',char(string(year(datetime(date)))),...
    '+-3monthsRAW.mat'],'Classifiers1yr','Classifiers2yr')
%% Bar plot
% Frequency of classifiers

% Find how many patient don't have class
f1yr=sum(~isnan(Classifiers1yr{:,2:end}));
f2yr=sum(~isnan(Classifiers2yr{:,2:end}));

figure('Name','Frequency of values per patient after surgery');
hold on
    g = vertcat(f1yr, f2yr);
    h = bar(g','EdgeColor','none','FaceColor','flat');
    h(1).FaceAlpha = 0.5;
    h(1).FaceColor = '#DD7A12';
    h(2).FaceColor = '#9800A0';
    h(2).FaceAlpha = 0.5;
    hline = refline([0 700]);
    hline.Color = 'w'; hline.LineWidth = 4;
    ax = gca;
    ax.FontSize = 14;
    xlabel('\fontsize{15} Outcome metrics',...
        'FontWeight','bold');
    ylabel('\fontsize{15} Frequency','FontWeight','bold');
    legend('\fontsize{12} 270-450 days','\fontsize{12} 630-810 days',...
        '\fontsize{12} threshold','Location','BestOutside');
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 28; %Width
    fig.Position(4) = 20; %Height
    xticklabels({'TWL','','WL'})
    print([subfold,'/',d,'/Plots/ClassifiersRAW'],'-deps')
    print([subfold,'/',d,'/Plots/ClassifiersRAWc'],'-depsc')
    fig.PaperOrientation = 'landscape';
    print([subfold,'/',d,'/Plots/ClassifiersRAW'], '-dpdf', '-fillpage')
hold off

clearvars hline g b;

close all
%%
% _Created by Aldo Arévalo_ 