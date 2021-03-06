%% Code to estimate Sucess or no-sucess class
% Sources: Pre-operativeValues2018 and FollowUpData

load([subfold,'/',d,'/PreOperativeValues.mat'])

load([subfold,'/',d,'/ImportedTables.mat'],'FollowUpData')

% Row's names
% patients_list = unique(GeneralData.PatientNr);

% Column's names
var_names = {'TWL','WL','EBMIL','EBWL','kg','BMI','m'};

%% Pre-allocate matrix
preallocate_nan = array2table(NaN(length(patients_list), length(var_names)+1));

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
    idx_patientPRE = [bmi_max.PatientCode] == patient_ID;
    BMI_patientPRE = bmi_max.Max(idx_patientPRE);
    
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
        
            
        % BMI
        idx_patientPRE = [height_mean.PatientCode] == patient_ID;
        height_patient = (height_mean.Mean(idx_patientPRE))/100;
        
            if isnan(weight_1yr) == 0
                BMI_1yr = weight_1yr/(height_patient^2);
            else
                BMI_1yr = NaN;
            end

            if isnan(weight_2yr) == 0
                BMI_2yr = weight_2yr/(height_patient^2);
            else
                 BMI_2yr = NaN;
            end
            
   Classifiers1yr.TWL(i) = (Weight_patientPRE - weight_1yr)/Weight_patientPRE;
   Classifiers1yr.WL(i) = Weight_patientPRE - weight_1yr;
   Classifiers1yr.EBMIL(i) = 100-(((BMI_1yr-25)/(BMI_patientPRE - 25))*100);
   Classifiers1yr.EBWL(i) = (Weight_patientPRE - weight_1yr)/(Weight_patientPRE-(25*(height_patient^2)));
   Classifiers1yr.kg(i) = weight_1yr;
   Classifiers1yr.BMI(i) = BMI_1yr;
   Classifiers1yr.m(i) = height_patient;
   
   Classifiers2yr.TWL(i) = (Weight_patientPRE - weight_2yr)/Weight_patientPRE;
   Classifiers2yr.WL(i) = Weight_patientPRE - weight_2yr;
   Classifiers2yr.EBMIL(i) = 100-(((BMI_2yr-25)/(BMI_patientPRE - 25))*100);
   Classifiers2yr.EBWL(i) = (Weight_patientPRE - weight_2yr)/(Weight_patientPRE-(25*(height_patient^2)));
   Classifiers2yr.kg(i) = weight_2yr;
   Classifiers2yr.BMI(i) = BMI_2yr;
   Classifiers2yr.m(i) = height_patient;
   
clearvars Weight_patientPRE weight_1yr weight_2yr patient_ID idx_patientFP ...
    idx_patientPRE height_patient BMI_patientPRE BMI_2yr BMI_1yr i weights_patients
end

% Delete negative values
%Classifiers2yr((Classifiers2yr.TWL<=0),:)=[];
%Classifiers1yr((Classifiers1yr.TWL<=0),:)=[];

% Convert negative values to zero
%varfun(@(x) x < 0, Classifiers1yr(:,{'TWL'}))
%varfun(@(x) x < 0, Classifiers2yr(:,{'TWL'}))

save([subfold,'/',d,'Classifiers2019+-3months.mat'],'Classifiers1yr','Classifiers2yr')
%% Bar plot
% Frequency of classifiers

% Find how many patient don't have class
f1yr=sum(~isnan(Classifiers1yr{:,2:end}));
f2yr=sum(~isnan(Classifiers2yr{:,2:end}));

figure('Name','Frequency of values per patient after surgery');
hold on
    g = vertcat(f1yr, f2yr);
    h = bar(g');
    h(1).FaceAlpha = 0.3;
    h(1).FaceColor = [0.1 0.1 1];
    h(2).FaceColor = [1 0.2 0.2];
    h(2).FaceAlpha = 0.3;
    hline = refline([0 700]);
    hline.Color = 'r';
    ax = gca;
    ax.FontSize = 12;
    xlabel('\fontsize{15} Classifiers. 1)%TWL, 2)WL, 3)%EBMIL, 4)%EBWL');
    ylabel('\fontsize{15} Frequency');
    legend('\fontsize{12} 270-450 days','\fontsize{12} 630-810 days','\fontsize{12} threshold');
hold off

clearvars hline g b;

%% Histograms
% Distribution of classifiers results
figure ('Name','Distribution of values after 1 yr and 2 yrs')
subplot(4,2,1)
    histogram(Classifiers1yr.TWL,'FaceAlpha',0.3);
    %xlim([0 Inf])
    ylabel('Frequency')
    title('Total Weight Loss %');
    legend('270-450 days');
        
subplot(4,2,2)
    histogram(Classifiers2yr.TWL,'FaceColor',[1 0.2 0.2]);
    %xlim([0 Inf])
    ylabel('Frequency')
    title('Total Weight Loss %');
    legend('630-810 days');

subplot(4,2,3)
    histogram(Classifiers1yr.WL,'FaceAlpha',0.3);
    %xlim([0 Inf])
    ylabel('Frequency')
    title('Absolute Weight Loss');
    xlabel('kg');
    legend('270-450 days');

subplot(4,2,4)
    histogram(Classifiers2yr.WL,'FaceColor',[1 0.2 0.2]);
    %xlim([0 Inf])
    ylabel('Frequency')
    title('Absolute Weight Loss');
    xlabel('kg');
    legend('630-810 days');

subplot(4,2,5)
    histogram(Classifiers1yr.EBMIL,'BinMethod', 'integer','FaceAlpha',0.3,...
        'EdgeColor','none');
    %xlim([0 Inf])
    ylabel('Frequency')
    title('Excess BMI Loss %');
    legend('270-450 days');

subplot(4,2,6)
    histogram(Classifiers2yr.EBMIL,'BinMethod', 'integer','EdgeColor','none',...
        'FaceColor',[1 0.2 0.2]);
    %xlim([0 Inf])
    ylabel('Frequency')
    title('Excess BMI Loss %');
    legend('630-810 days');
    
subplot(4,2,7)
    histogram(Classifiers1yr.EBWL,'FaceAlpha',0.3);
    %xlim([0 Inf])
    ylabel('Frequency')
    title('Excess Body Weight Loss %');
    legend('270-450 days');

subplot(4,2,8)
    histogram(Classifiers2yr.EBWL,'FaceAlpha',0.3);
    %xlim([0 Inf])
    ylabel('Frequency')
    title('Excess Body Weight Loss %');
    legend('630-810 days');
    
fig = gcf;
fig.Units = 'centimeter';
fig.Position(1)= 5; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 3; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 35; %Width
fig.Position(4) = 22; %Height

print([subfold,'/',d,'/Plots/outputs'],'-deps')
print([subfold,'/',d,'/Plots/outputsc'],'-depsc')
print([subfold,'/',d,'/Plots/outputs'], '-dpdf','-bestfit')

close all

%%
% _Created by Aldo Ar�valo_ 
