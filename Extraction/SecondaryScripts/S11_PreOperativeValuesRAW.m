%% Code to estimate WEIGHT, BMI, HEIGHT and WAIST CIRCUMFERENCE
% Sources: DATOScreening, Physical Measures and Pre-Surgical Screening

load([subfold,'/',d,'/RawData/PhMsDataRAW.mat'])

load([subfold,'/',d,'/RawData/PPOSdataRAW2.mat'])

load([subfold,'/',d,'/RawData/ScreeningDataRAW.mat'])

load([subfold,'/',d,'/RawData/FollowUpDataRAW.mat'], 'FollowUpData')

load([subfold,'/',d,'/RawData/LabDataRAW.mat'])

load([subfold,'/',d,'/RawData/GeneralDataRAW.mat'])
total_patientsID = unique(GeneralData.PatientNr);

% Row's names
patients_list = total_patientsID;

%% Create LabResultsPerPatient for RAW data

col_2_add = array2table(NaN(height(LabResults),3));
col_2_add.Properties.VariableNames = {'RefLow','RefUp','StandardizedResult'};

LabResults = [LabResults col_2_add];

% Split the field 'refwaarde' to two seperate fields for lower and upper
    % boundaries. This column contains the normal ranges for an specific
    % biomarker/lab test according to age, sex and gender of each patient

for j=1:length(total_patientsID)
    % Get patient ID for which data will be collected
    patientID = total_patientsID(j);
    % Get index vector where the data for this patient is located in
    % LabResults table
    patientidx = LabResults.PatientCode == patientID;
    % Gather data for the j patient
    A = LabResults(patientidx,:);
    % Split string of each row where '-' is located
    split = cellfun(@(x) regexp(char(x),'-', 'split'),...
        table2array(LabResults(patientidx,'refwaarde')),...
        'UniformOutput',false);
    % Create vector to pre-allocate values
    ref = NaN(height(LabResults(patientidx,:)),2);
    
    % Find cases where the following constraints are fullfilled
    for i=1:height(LabResults(patientidx,:))
        
        if ~cellfun('isempty',(strfind(split{i,1},'>')))
            % Split one-sided (with lower) reference intervals
            splitnew = strrep(split{i,1},'>','');
            ref(i,1) = str2double(splitnew);
            
        elseif ~cellfun('isempty',(strfind(split{i,1},'<')))
            % Split one-sided (with upper) reference intervals
            splitnew = strrep(split{i,1},'<','');
            ref(i,1) = 0;
            ref(i,2) = str2double(splitnew);
            
        elseif sum(size(split{i,1}))> 2
            % Split double-sided reference intervals
            ref(i,1) = str2double(split{i,1}{1,1});
            ref(i,2) = str2double(split{i,1}{1,2});
        end
        
        clearvars splitnew
        
        % Check the constraints for the following lab tests
        % Fill in reference values for lab results that are manually computed
        if and(isnan(ref(i,1)),isnan(ref(i,2)))
            if strcmp(A{i,'bepcode'},'BCH014') % Cholesterol/HDL ratio
                ref(i,1) = 0;
                ref(i,2) = 8;
            elseif strcmp(A{i,'bepcode'},'BTR003') % Triglycerides (mM)
                ref(i,1) = 0;
                ref(i,2) = 1.8;
            elseif strcmp(A{i,'bepcode'},'BCH004') % Cholesterol (mM)
                ref(i,1) = 0;
                ref(i,2) = 5;
            elseif strcmp(A{i,'bepcode'},'BHD000') % HDL - Cholesterol (mM)
                ref(i,1) = 1;
            elseif strcmp(A{i,'bepcode'},'BLD002') % LDL - Cholesterol (mM)
                    ref(i,1) = 0;
                    ref(i,2) = 2.5;
            elseif strcmp(A{i,'bepcode'},'BAL015') % ALAT (U/L)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    ref(i,1) = 0;
                    ref(i,2) = 45;
                elseif GeneralData{t,'Sex'} == 0 % Females
                    ref(i,1) = 0;
                    ref(i,2) = 34;
                end
                clearvars t
            elseif strcmp(A{i,'bepcode'},'BAS002') % ASAT (U/L)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    ref(i,1) = 0;
                    ref(i,2) = 35;
                elseif GeneralData{t,'Sex'} == 0 % Females
                    ref(i,1) = 0;
                    ref(i,2) = 31;
                end
                clearvars t
            elseif strcmp(A{i,'bepcode'},'BCR002') % CRP (mg/L)
                ref(i,1) = 0;
                ref(i,2) = 6;
            elseif strcmp(A{i,'bepcode'},'BCO017') % Cortisol (nM)
                ref(i,1) = 130;
                ref(i,2) = 540;
            elseif strcmp(A{i,'bepcode'},'BGT001') % gamma GT (U/L)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    ref(i,1) = 0;
                    ref(i,2) = 55;
                elseif GeneralData{t,'Sex'} == 0 % Females
                    ref(i,1) = 0;
                    ref(i,2) = 38;
                end
                clear t
            elseif strcmp(A{i,'bepcode'},'BGL003') % Glucose (blood) (mM)
                ref(i,1) = 0;
                ref(i,2) = 7.8;
            elseif strcmp(A{i,'bepcode'},'BIN000') % Insulin (mU/L)
                ref(i,1) = 0;
                ref(i,2) = 17;
            elseif strcmp(A{i,'bepcode'},'BIN008') % INR (PT)
                ref(i,1) = 2.0;
                ref(i,2) = 3.5;
            elseif strcmp(A{i,'bepcode'},'BIN001') % Inhalant allergy screening
                % Boolean normal = negative (0)
                ref(i,1) = 0;
                ref(i,2) = 1;
                
            elseif strcmp(A{i,'bepcode'},'UBI000') % Bilirubin in urine
                % Boolean normal = negative (0)
                ref(i,1) = 0;
                ref(i,2) = 1;
                
            elseif strcmp(A{i,'bepcode'},'BFE006') % Ferritine (ug/L)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    ref(i,1) = 30;
                    ref(i,2) = 400;
                elseif GeneralData{t,'Sex'} == 0 % Females
                    ref(i,1) = 13;
                    ref(i,2) = 200;
                end
                clear t
            elseif strcmp(A{i,'bepcode'},'BHC004') % hCG (+Beta-subunits) (U/L)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    ref(i,1) = 0;
                    ref(i,2) = 2.6;
                elseif GeneralData{t,'Sex'} == 0 % Females
                    ref(i,1) = 0;
                    ref(i,2) = 5.3;
                end
                clear t
            elseif strcmp(A{i,'bepcode'},'BAN037') % ANA anti-nuclear antibodies
                % Boolean normal = negative (0)
                ref(i,1) = 0;
                ref(i,2) = 1;
            elseif strcmp(A{i,'bepcode'},'BIJ001') % Iron (uM)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    ref(i,1) = 14;
                    ref(i,2) = 35;
                elseif GeneralData{t,'Sex'} == 0 % Females
                    ref(i,1) = 10;
                    ref(i,2) = 25;
                end
                clear t
            elseif strcmp(A{i,'bepcode'},'UKR002') % Creatinine in urine (mM)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    ref(i,1) = 3;
                    ref(i,2) = 29;
                elseif GeneralData{t,'Sex'} == 0 % Females
                    ref(i,1) = 2;
                    ref(i,2) = 28;
                end
                clear t
            elseif strcmp(A{i,'bepcode'},'BAM012') % Amylase (U/L)
                ref(i,1) = 0;
                ref(i,2) = 100;
            elseif strcmp(A{i,'bepcode'},'BBE001') % Sedimentation of erythrocytes
                ref(i,1) = 0;
                ref(i,2) = 6;
            elseif strcmp(A{i,'bepcode'},'BBI000') % Bilirubin conjugated (uM)
                ref(i,1) = 0;
                ref(i,2) = 5;
            elseif strcmp(A{i,'bepcode'},'UEI001') % Protein in urine
                % Boolean normal = negative (0)
                ref(i,1) = 0;
                ref(i,2) = 1;
            elseif strcmp(A{i,'bepcode'},'UGL000') % Glucose in urine
                % Boolean normal = negative (0)
                ref(i,1) = 0;
                ref(i,2) = 1;
            elseif strcmp(A{i,'bepcode'},'UKE000') % Ketonic compounds
                % Boolean normal = negative (0)
                ref(i,1) = 0;
                ref(i,2) = 1;
            elseif strcmp(A{i,'bepcode'},'UNI000') % Nitrite in urine
                % Boolean normal = negative (0)
                ref(i,1) = 0;
                ref(i,2) = 1;
            elseif strcmp(A{i,'bepcode'},'BBA002') % Basophyles (cells/nL)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    if GeneralData{t,'Age'} >= 50.0
                        ref(i,1) = 0;
                        ref(i,2) = 20;
                    else % < 50 years old
                        ref(i,1) = 0;
                        ref(i,2) = 15;    
                    end
                end
                if GeneralData{t,'Sex'} == 0 % Females
                   if GeneralData{t,'Age'} >= 50.0
                        ref(i,1) = 0;
                        ref(i,2) = 30;
                   else % < 50 years old
                        ref(i,1) = 0;
                        ref(i,2) = 20;    
                   end
                end
                clear t
            end
        end
        
        if and(isnan(ref(i,1)),~isnan(ref(i,2)))
            if strcmp(A{i,'bepcode'},'BAM012') % Amylase (U/L)
                ref(i,1) = 0;
                ref(i,2) = 100;
            elseif strcmp(A{i,'bepcode'},'BBA002') % Basophyles (cells/nL)
                ref(i,1) = 0;
                ref(i,2) = 0.2;
            elseif strcmp(A{i,'bepcode'},'BCR002') % CRP (mg/L)
                ref(i,1) = 0;
                ref(i,2) = 6;
            elseif strcmp(A{i,'bepcode'},'BBE001') % Sedimentation of erythrocytes
                ref(i,1) = 0;
                ref(i,2) = 6;
            elseif strcmp(A{i,'bepcode'},'BBI000') % Bilirubin conjugated (uM)
                ref(i,1) = 0;
                ref(i,2) = 5;
            elseif strcmp(A{i,'bepcode'},'OBA000') % Base excess in arterial blood(mM)
                ref(i,1) = -2;
                ref(i,2) = 3;
            elseif strcmp(A{i,'bepcode'},'UBI000') % Bilirubin in urine
                % Boolean normal = negative (0)
                ref(i,1) = 0;
                ref(i,2) = 1;
            elseif strcmp(A{i,'bepcode'},'UEI001') % Protein in urine
                % Boolean normal = negative (0)
                ref(i,1) = 0;
                ref(i,2) = 1;
            elseif strcmp(A{i,'bepcode'},'UKR002') % Creatinine in urine (mM)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    ref(i,1) = 3;
                    ref(i,2) = 29;
                elseif GeneralData{t,'Sex'} == 0 % Females
                    ref(i,1) = 2;
                    ref(i,2) = 28;
                end
                clear t    
            elseif strcmp(A{i,'bepcode'},'BBA002') % Basophyles (cells/nL)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    if GeneralData{t,'Age'} >= 50.0
                        ref(i,1) = 0;
                        ref(i,2) = 20;
                    else % < 50 years old
                        ref(i,1) = 0;
                        ref(i,2) = 15;    
                    end
                end
                if GeneralData{t,'Sex'} == 0 % Females
                   if GeneralData{t,'Age'} >= 50.0
                        ref(i,1) = 0;
                        ref(i,2) = 30;
                   else % < 50 years old
                        ref(i,1) = 0;
                        ref(i,2) = 20;  
                   end
                end
                clear t
            elseif strcmp(A{i,'bepcode'},'BHC004') % hCG (+Beta-subunits) (U/L)
                t = find(A{i,'PatientCode'} == table2array(GeneralData(:,'PatientNr')));
                if GeneralData{t,'Sex'} == 1 % Males
                    ref(i,1) = 0;
                    ref(i,2) = 2.6;
                elseif GeneralData{t,'Sex'} == 0 % Females
                    ref(i,1) = 0;
                    ref(i,2) = 5.3;
                end
                clear t
            end
        end
    end

ReferenceValues = array2table(ref);
ReferenceValues.Properties.VariableNames = {'RefLow','RefUp'};

% Allocate in original dataset
LabResults(patientidx,{'RefLow'}) = ReferenceValues(:,{'RefLow'});
LabResults(patientidx,{'RefUp'}) = ReferenceValues(:,{'RefUp'});

clear split ref i col_2_add
end

% Extract bepcodes in a cell array
var_names={'BHE001','BHE000','BER002','BMC000','BMC002',...
    'BTR006','BLE001','BGL003','BBI001','BAS002','BAL015','BLD004','BAL014',...
    'BGT001','BUR002','BKR000b','BKR015','BKA000','BNA001','BCA006b',...
    'BFO002','BAL002','BCR002','BCH004','BHD000','BCH014','BLD002',...
    'BPR012','BHE030','BIN000','BCP002','BPA016','BTS003','BCO017',...
    'BIJ001','BFE006','BFO000','BMA002'};

data_struct.Result = {};
data_struct.days_from_ok = NaN;

ColumSruct.PatientCode = NaN;
for i = 1:length(var_names)
    ColumSruct.(var_names{i}) = data_struct;
end

list_PatientID = unique(LabResults.PatientCode);

%% Pre-allocate structs for all patients
LabResultsPerPatient = repmat(ColumSruct,length(list_PatientID),1);
clearvars i;

for i=1:length(list_PatientID)
    patient_ID = list_PatientID(i);
    LabResultsPerPatient(i).PatientCode = patient_ID;
    idx_patient_ID = LabResults.PatientCode == patient_ID;
    s = LabResults(idx_patient_ID,:);
    var_names = unique(s.bepcode);
    
    for j = 1:length(var_names)
        var_name = var_names{j};
        idx_var = strcmp(s.bepcode,var_name);
        
        vars = table2array(s(idx_var,{'uitslagNUM'}));
        n = size(vars,1);
        
        analysis_list = {};
        for k = 1:n
            analysis_list{k,1} = vars(k,~isnan(vars(k,:)));
        end
        days_from_ok_list = s.days_from_ok(idx_var);
        
        analysisList_empty = find(cellfun(@(x) isempty(x),analysis_list));
        analysis_list(analysisList_empty,:) = [];
        days_from_ok_list(analysisList_empty,:) = [];
        
            if isempty(analysis_list) == 0
                LabResultsPerPatient(i).(var_name).Result = cell2mat(analysis_list);
                LabResultsPerPatient(i).(var_name).days_from_ok = days_from_ok_list;       
            else
                LabResultsPerPatient(i).(var_name).Result = NaN;
                LabResultsPerPatient(i).(var_name).days_from_ok = NaN;
            end
    end
    clear idx_var vars n s var_names idx_patient_ID var_name ...
         days_from_ok_list k j i data_struct analysis_list patient_ID
end

clear data_struct 

%% Obtain screening date per patient RAW data

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
                        LabPackTimesPatient.(var_name)(i) = median(unique_days(vidx),'omitnan');
                    else
                        LabPackTimesPatient.(var_name)(i) = NaN;
                    end
                otherwise % >= 3
                    vidx = arrayfun(@(x) x >= -365, unique_days);
                    if sum(vidx) >= 1
                        LabPackTimesPatient.(var_name)(i) = median(unique_days(vidx),'omitnan');
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

clearvars A
%%  Column's names
% BRI = Body Roundness Index
% ABSI = A Body Shape Index
% TBFM = Total Body Fat Mass

var_names = {'Weight','BMI','Height','WaistCircum','BRI','ABSI','TBFM'};

%% Fitting function
%ALLFITDIST Fit all valid parametric probability distributions to data.
%   [D PD] = ALLFITDIST(DATA) fits all valid parametric probability
%   distributions to the data in vector DATA, and returns a struct D of
%   fitted distributions and parameters and a struct of objects PD
%   representing the fitted distributions. PD is an object in a class
%   derived from the ProbDist class.
%
%   [...] = ALLFITDIST(DATA,SORTBY) returns the struct of valid distributions
%   sorted by the parameter SORTBY
%        NLogL - Negative of the log likelihood
%        BIC - Bayesian information criterion (default)
%        AIC - Akaike information criterion
%        AICc - AIC with a correction for finite sample sizes

% In statistics, the Bayesian information criterion (BIC) or Schwarz 
% criterion (also SBC, SBIC) is a criterion for model selection among a 
% finite set of models; the model with the lowest BIC is preferred. It is
% based, in part, on the likelihood function and it is closely related to
% the Akaike information criterion (AIC).

%Maximum number of distributions to include
max_num_dist=Inf;  %All valid distributions

%% Pre-allocate matrix
prealloc_num = array2table(NaN(length(patients_list), length(var_names)+1));

PreOpMean = prealloc_num;
PreOpMin = prealloc_num;
PreOpMax = prealloc_num;
PreOpNM = prealloc_num;
PreOpLast = prealloc_num;
PreOpMedian = prealloc_num;
PreOpNMPrior = prealloc_num;
PreOpNMTotal = prealloc_num;

% Add column names
col_names = horzcat('PatientCode',var_names);

PreOpMean.Properties.VariableNames = col_names;
PreOpMin.Properties.VariableNames = col_names;
PreOpMax.Properties.VariableNames = col_names;
PreOpNM.Properties.VariableNames = col_names;
PreOpLast.Properties.VariableNames = col_names;
PreOpMedian.Properties.VariableNames = col_names;
PreOpNMPrior.Properties.VariableNames = col_names;
PreOpNMTotal.Properties.VariableNames = col_names;

% For estimating output variable
weight_max = array2table(NaN(length(patients_list), 2));
weight_max.Properties.VariableNames = {'PatientCode','Max'};

bmi_max = array2table(NaN(length(patients_list), 2));
bmi_max.Properties.VariableNames = {'PatientCode','Max'};

height_mean = array2table(NaN(length(patients_list), 2));
height_mean.Properties.VariableNames = {'PatientCode','Mean'};

% Indexes of variables
    % Weight
    idx_weight_PPOS = strcmp(PPOSdata.questionID, 'vID148748');
        aux1 = [PhMsData.vraagID] == 111972;
        aux2 = [PhMsData.vraagID] == 113882;
        aux3 = [PhMsData.vraagID] == 118486;
        aux4 = [PhMsData.vraagID] == 121101;
        aux5 = [PhMsData.vraagID] == 121008;
    idx_weight_PhMs = logical(sum(horzcat(aux1, aux2, aux3, aux4,aux5),2));
    
    clearvars aux1 aux2 aux3 aux4 aux5 total_patientsID
    
    % Check distribution
    Ant = table2cell(PhMsData(idx_weight_PhMs,{'antwoord'}));
    Ant = cellfun(@str2double,Ant);
    WeightALL=vertcat(table2array(PPOSdata(idx_weight_PPOS,{'Value01'})),...
        Ant,table2array(PhMsData(idx_weight_PhMs,{'value01'})));
    WeightALL(isnan(WeightALL))=[];
    WeightALL(~any(WeightALL,2))=[];
    [D,PD] = allfitdist(WeightALL);
    
    mWeight = NaN(3,1);
    sWeight = NaN(3,1);
    DWeight = D;
    
    for i=1:3
        s = std(PD{i}); sWeight(i,1) = s;
        m = mean(PD{i}); mWeight(i,1) = m;
        figure('Name',['Weight ',D(i).DistName])
        nbins = max(min(length(WeightALL)./10,100),50);
        xi = linspace(min(WeightALL),max(WeightALL),nbins);
        dx = mean(diff(xi));
        xi2 = linspace(min(WeightALL),max(WeightALL),nbins*10)';
        fi = histc(WeightALL,xi-dx);
        fi = fi./sum(fi)./dx;
        ys = cellfun(@(PD) pdf(PD,xi2),PD(i),'UniformOutput',0);
        ys = cat(2,ys{:});
        b1 = bar(xi,fi,'FaceColor',[160 188 254]/255,'EdgeColor','k');
        hold on
        p1 = plot(xi2,ys,'LineWidth',1.5);
        ax = gca;
        xL=get(ax,'Ylim');
        h1=line([m m], xL, 'Color','r','LineStyle','--','LineWidth',1); % mean
        h2=line([m+3*s m+3*s], xL, 'Color','m','LineStyle','--','LineWidth',1,'MarkerEdgeColor','k'); % +ST dev
        h3=line([m-3*s m-3*s], xL, 'Color','m','LineStyle','--','LineWidth',1,'MarkerEdgeColor','k'); % +ST dev
        h4=line([m+4*s m+4*s], xL, 'Color','b','LineStyle','--','LineWidth',1); % +ST dev
        h5=line([m-4*s m-4*s], xL, 'Color','b','LineStyle','--','LineWidth',1); % +ST dev
        legend([b1 p1 h1 h3 h5],'Real data','PDF','mean','3SD','4SD','Location','bestoutside')
        title(D(i).DistName)
        xlabel('Weight (kg)','FontSize',18,'FontWeight','bold') % x-axis label
        ylabel('Probability Density','FontSize',18,'FontWeight','bold') % y-axis label
        ax.FontSize = 14; ax.FontWeight = 'bold';
        hold off
        
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1)= 5; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 3; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 35; %Width
        fig.Position(4) = 22; %Height
        fig.PaperOrientation = 'landscape';
        
%         print([subfold,'/',d,'/Plots/Weight',D(i).DistName],'-deps')
%         print([subfold,'/',d,'/Plots/Weight',D(i).DistName,'c'],'-depsc')
%         print([subfold,'/',d,'/Plots/Weight',D(i).DistName], '-dpdf','-bestfit')
    end
    
    close all
    
    % BMI
    idx_BMI_PPOS = strcmp(PPOSdata.questionID, 'vID148749');
        aux1 = [PhMsData.vraagID] == 116086;
        aux2 = [PhMsData.vraagID] == 118678;
    idx_BMI_PhMs = logical(sum(horzcat(aux1, aux2),2));
    
    clearvars aux1 aux2
    
    % Check distribution
    Ant = table2cell(PPOSdata(idx_BMI_PPOS,{'antwoord'}));
    Ant = cellfun(@str2double,Ant);
    Ant2 = table2cell(PhMsData(idx_BMI_PhMs,{'antwoord'}));
    Ant2 = cellfun(@str2double,Ant2);
                
    BMIall = vertcat(Ant,Ant2);
    BMIall(isnan(BMIall))=[];
    BMIall(~any(BMIall,2))=[];
    [D,PD] = allfitdist(BMIall);
    
    mBMI = NaN(3,1);
    sBMI = NaN(3,1);
    DBMI = D;
    
    for i=1:3
        s = std(PD{i}); sBMI(i,1) = s;
        m = mean(PD{i}); mBMI(i,1) = m;
        figure('Name',['BMI ',D(i).DistName])
        nbins = max(min(length(BMIall)./10,100),50);
        xi = linspace(min(BMIall),max(BMIall),nbins);
        dx = mean(diff(xi));
        xi2 = linspace(min(BMIall),max(BMIall),nbins*10)';
        fi = histc(BMIall,xi-dx);
        fi = fi./sum(fi)./dx;
        ys = cellfun(@(PD) pdf(PD,xi2),PD(i),'UniformOutput',0);
        ys = cat(2,ys{:});
        b1 = bar(xi,fi,'FaceColor',[160 188 254]/255,'EdgeColor','k');
        hold on
        p1 = plot(xi2,ys,'LineWidth',1.5);
        ax = gca;
        xL=get(ax,'Ylim');
        h1=line([m m], xL, 'Color','r','LineStyle','--','LineWidth',1); % mean
        h2=line([m+3*s m+3*s], xL, 'Color','m','LineStyle','--','LineWidth',1,'MarkerEdgeColor','k'); % +ST dev
        h3=line([m-3*s m-3*s], xL, 'Color','m','LineStyle','--','LineWidth',1,'MarkerEdgeColor','k'); % +ST dev
        h4=line([m+4*s m+4*s], xL, 'Color','b','LineStyle','--','LineWidth',1); % +ST dev
        h5=line([m-4*s m-4*s], xL, 'Color','b','LineStyle','--','LineWidth',1); % +ST dev
        legend([b1 p1 h1 h3 h5],'Real data','PDF','mean','3SD','4SD','Location','bestoutside')
        title(D(i).DistName)
        xlabel('BMI (kg/m^{2})','FontSize',18,'FontWeight','bold') % x-axis label
        ylabel('Probability Density','FontSize',18,'FontWeight','bold') % y-axis label
        ax.FontSize = 14; ax.FontWeight = 'bold';
        hold off
        
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1)= 5; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 3; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 35; %Width
        fig.Position(4) = 22; %Height
        fig.PaperOrientation = 'landscape';
        
%         print([subfold,'/',d,'/Plots/BMI',D(i).DistName],'-deps')
%         print([subfold,'/',d,'/Plots/BMI',D(i).DistName,'c'],'-depsc')
%         print([subfold,'/',d,'/Plots/BMI',D(i).DistName], '-dpdf','-bestfit')
    end
    
    close all
    
    % Height
    idx_height_PPOS = strcmp(PPOSdata.questionID, 'vID148747');
    idx_height_PhMs=table2array(varfun(@(x) (x==116084 | x==118487 | ...
        x== 121100 | x== 132959),PhMsData(:,{'vraagID'})));
    
    % Convert to cm ScreeningData.lengte colum
    [ScreeningData.lengte]= ScreeningData.lengte*100;
    
    % Check distribution
    Ant = table2cell(PhMsData(idx_height_PhMs,{'antwoord'}));
    Ant = cellfun(@str2double,Ant);
                
    HeightALL = vertcat(table2array(PPOSdata(idx_height_PPOS,{'Value01'})),...
        Ant,table2array(PhMsData(idx_height_PhMs,{'value01'})),...
        ScreeningData.lengte);
    HeightALL(isnan(HeightALL))=[];
    HeightALL(HeightALL<58)=[];
    [D,PD] = allfitdist(HeightALL);
    
    mHeight = NaN(3,1);
    sHeight = NaN(3,1);
    DHeight = D;

    for i=1:3
        s = std(PD{i}); sHeight(i,1) = s;
        m = mean(PD{i}); mHeight(i,1) = m;
        figure('Name',['Height ',D(i).DistName])
        nbins = max(min(length(HeightALL)./10,100),50);
        xi = linspace(min(HeightALL),max(HeightALL),nbins);
        dx = mean(diff(xi));
        xi2 = linspace(min(HeightALL),max(HeightALL),nbins*10)';
        fi = histc(HeightALL,xi-dx);
        fi = fi./sum(fi)./dx;
        ys = cellfun(@(PD) pdf(PD,xi2),PD(i),'UniformOutput',0);
        ys = cat(2,ys{:});
        b1 = bar(xi,fi,'FaceColor',[160 188 254]/255,'EdgeColor','k');
        hold on
        p1 = plot(xi2,ys,'LineWidth',1.5);
        ax = gca;
        xL=get(ax,'Ylim');
        h1=line([m m], xL, 'Color','r','LineStyle','--','LineWidth',1); % mean
        h2=line([m+3*s m+3*s], xL, 'Color','m','LineStyle','--','LineWidth',1,'MarkerEdgeColor','k'); % +ST dev
        h3=line([m-3*s m-3*s], xL, 'Color','m','LineStyle','--','LineWidth',1,'MarkerEdgeColor','k'); % +ST dev
        h4=line([m+4*s m+4*s], xL, 'Color','b','LineStyle','--','LineWidth',1); % +ST dev
        h5=line([m-4*s m-4*s], xL, 'Color','b','LineStyle','--','LineWidth',1); % +ST dev
        legend([b1 p1 h1 h3 h5],'Real data','PDF','mean','3SD','4SD','Location','bestoutside')
        title(D(i).DistName)
        xlabel('Height (cm)','FontSize',18,'FontWeight','bold') % x-axis label
        ylabel('Probability Density','FontSize',18,'FontWeight','bold') % y-axis label
        ax.FontSize = 14; ax.FontWeight = 'bold';
        hold off
        
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1)= 5; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 3; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 35; %Width
        fig.Position(4) = 22; %Height
        fig.PaperOrientation = 'landscape';
        
%         print([subfold,'/',d,'/Plots/Height',D(i).DistName],'-deps')
%         print([subfold,'/',d,'/Plots/Height',D(i).DistName,'c'],'-depsc')
%         print([subfold,'/',d,'/Plots/Height',D(i).DistName], '-dpdf','-bestfit')
    end
    
    close all
    
    % Waist circumference
    idx_buik_PhMs = [PhMsData.vraagID] == 118489;
    WCall = table2array(PhMsData(idx_buik_PhMs,{'value01'}));
    WCall(isnan(WCall))=[];
    WCall(~any(WCall,2))=[];
    
    [D,PD] = allfitdist(WCall);
    
    mWC = NaN(3,1);
    sWC = NaN(3,1);
    DWC = D;

    for i=1:3
        s = std(PD{i}); sWC(i,1) = s;
        m = mean(PD{i}); mWC(i,1) = m;
        figure('Name',['Height ',D(i).DistName])
        nbins = max(min(length(WCall)./10,100),50);
        xi = linspace(min(WCall),max(WCall),nbins);
        dx = mean(diff(xi));
        xi2 = linspace(min(WCall),max(WCall),nbins*10)';
        fi = histc(WCall,xi-dx);
        fi = fi./sum(fi)./dx;
        ys = cellfun(@(PD) pdf(PD,xi2),PD(i),'UniformOutput',0);
        ys = cat(2,ys{:});
        b1 = bar(xi,fi,'FaceColor',[160 188 254]/255,'EdgeColor','k');
        hold on
        p1 = plot(xi2,ys,'LineWidth',1.5);
        ax = gca;
        xL=get(ax,'Ylim');
        h1=line([m m], xL, 'Color','r','LineStyle','--','LineWidth',1); % mean
        h2=line([m+3*s m+3*s], xL, 'Color','m','LineStyle','--','LineWidth',1,'MarkerEdgeColor','k'); % +ST dev
        h3=line([m-3*s m-3*s], xL, 'Color','m','LineStyle','--','LineWidth',1,'MarkerEdgeColor','k'); % +ST dev
        h4=line([m+4*s m+4*s], xL, 'Color','b','LineStyle','--','LineWidth',1); % +ST dev
        h5=line([m-4*s m-4*s], xL, 'Color','b','LineStyle','--','LineWidth',1); % +ST dev
        legend([b1 p1 h1 h3 h5],'Real data','PDF','mean','3SD','4SD','Location','bestoutside')
        title(D(i).DistName)
        xlabel('Waist Circumference (cm)','FontSize',18,'FontWeight','bold') % x-axis label
        ylabel('Probability Density','FontSize',18,'FontWeight','bold') % y-axis label
        ax.FontSize = 14; ax.FontWeight = 'bold';
        hold off
        
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1)= 5; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 3; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 35; %Width
        fig.Position(4) = 22; %Height
        fig.PaperOrientation = 'landscape';
        
%         print([subfold,'/',d,'/Plots/WC',D(i).DistName],'-deps')
%         print([subfold,'/',d,'/Plots/WC',D(i).DistName,'c'],'-depsc')
%         print([subfold,'/',d,'/Plots/WC',D(i).DistName], '-dpdf','-bestfit')
    end
    
    close all
    
clearvars ax b1 D dx fi fig h1 h2 h3 h4 h5 i m max_num_dist nbins p1 PD ...
    s xi xi2 xL ys

%% Delete unadmissible values over 4SD
    % 41.3817 kg < Weight > 213.7769 kg
    Ant = varfun(@str2double,PhMsData(:,{'antwoord'}));
    PhMsData(:,{'antwoord'}) = [];
    PhMsData(:,{'antwoord'}) = Ant;
    
    g1 = table2array(varfun(@(x) (x<(mWeight(1)-4*sWeight(1)) | x>(mWeight(1)+4*sWeight(1))),PhMsData(:,{'antwoord', 'value01', 'value02'})));
    rows_antwoord_w = idx_weight_PhMs & g1(:,1);
    rows_value01_w = idx_weight_PhMs & g1(:,2);
    rows_value02_w = idx_weight_PhMs & g1(:,3);
    PhMsData.antwoord(rows_antwoord_w)=NaN;
    PhMsData.value01(rows_value01_w)=NaN;
    PhMsData.value02(rows_value02_w)=NaN;
    
    h1 = table2array(varfun(@(x) (x<(mWeight(1)-4*sWeight(1)) | x>(mWeight(1)+4*sWeight(1))),PPOSdata(:,{'antwoord', 'Value01'})));
    rows_antwoord2_w = idx_weight_PPOS & h1(:,1);
    rows_Value01_w = idx_weight_PPOS & h1(:,2);
    PPOSdata.antwoord(rows_antwoord2_w)=NaN;
    PPOSdata.Value01(rows_Value01_w)=NaN;
    
    hooggew_70 = ScreeningData.hooggew<(mWeight(1)-4*sWeight(1));
    ScreeningData.hooggew(hooggew_70)= NaN;
    hooggew_220 = ScreeningData.hooggew>(mWeight(1)+4*sWeight(1));
    ScreeningData.hooggew(hooggew_220)= NaN;
    
    gewicht_70 = ScreeningData.gewicht<(mWeight(1)-4*sWeight(1));
    ScreeningData.gewicht(gewicht_70)= NaN;
    gewicht_220 = ScreeningData.gewicht>(mWeight(1)+4*sWeight(1));
    ScreeningData.gewicht(gewicht_220)= NaN;
    
    % Check ditribution
    WeightALL=vertcat(table2array(PPOSdata(idx_weight_PPOS,{'Value01'})),...
        table2array(PhMsData(idx_weight_PhMs,{'antwoord'})),...
        table2array(PhMsData(idx_weight_PhMs,{'value01'})));
    WeightALL(isnan(WeightALL))=[];
    WeightALL(~any(WeightALL,2))=[];
    
    figure('Name','Weight final')
    b1=histogram(WeightALL,'DisplayName','Weight (kg)','FaceColor','#EDB120',...
            'EdgeColor','none');
    hold on
    ax = gca;
    xL=get(ax,'Ylim');
    h1=line([mWeight(1) mWeight(1)], xL, 'Color','r','LineStyle','--','LineWidth',2); % mean
    h2=line([mWeight(1)+3*sWeight(1) mWeight(1)+3*sWeight(1)], xL, 'Color','#8900BF','LineStyle','--','LineWidth',2,'MarkerEdgeColor','k'); % +ST dev
    h3=line([mWeight(1)-3*sWeight(1) mWeight(1)-3*sWeight(1)], xL, 'Color','#8900BF','LineStyle','--','LineWidth',2,'MarkerEdgeColor','k'); % +ST dev
    h4=line([mWeight(1)+4*sWeight(1) mWeight(1)+4*sWeight(1)], xL, 'Color','b','LineStyle','--','LineWidth',2); % +ST dev
    h5=line([mWeight(1)-4*sWeight(1) mWeight(1)-4*sWeight(1)], xL, 'Color','b','LineStyle','--','LineWidth',2); % +ST dev
    legend([b1 h1 h3 h5],'Real data','mean','3SD','4SD','Location','bestoutside')
    title('Weight final distribution')
    xlabel('Weight (kg)','FontSize',18,'FontWeight','bold') % x-axis label
    ylabel('Frequency','FontSize',18,'FontWeight','bold') % y-axis label
    ax.FontSize = 14; ax.FontWeight = 'bold';
    hold off
        
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1)= 5; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 3; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 35; %Width
    fig.Position(4) = 22; %Height
    fig.PaperOrientation = 'landscape';
        
%     print([subfold,'/',d,'/Plots/Weightfinal'],'-deps')
%     print([subfold,'/',d,'/Plots/Weightfinal','c'],'-depsc')
%     print([subfold,'/',d,'/Plots/Weightfinal'], '-dpdf','-bestfit')
    
    % Height
    g2 = table2array(varfun(@(x) (x<(mHeight(1)-4*sHeight(1))),PhMsData(:,{'antwoord', 'value01', 'value02'})));
    rows_antwoord_h = idx_height_PhMs & g2(:,1);
    rows_value01_h = idx_height_PhMs & g2(:,2);
    rows_value02_h = idx_height_PhMs & g2(:,3);
    PhMsData.antwoord(rows_antwoord_h)=NaN;
    PhMsData.value01(rows_value01_h)=NaN;
    PhMsData.value02(rows_value02_h)=NaN;
    
    h2 = table2array(varfun(@(x) x<(mHeight(1)-4*sHeight(1)),PPOSdata(:,{'antwoord', 'Value01'})));
    rows_antwoord2_h = idx_height_PPOS & h2(:,1);
    rows_Value01_h = idx_height_PPOS & h2(:,2);
    PPOSdata.antwoord(rows_antwoord2_h)=NaN;
    PPOSdata.Value01(rows_Value01_h)=NaN;
    
    % Check distribution
    HeightALL = vertcat(table2array(PPOSdata(idx_height_PPOS,{'Value01'})),...
        table2array(PhMsData(idx_height_PhMs,{'antwoord'})),...
        table2array(PhMsData(idx_height_PhMs,{'value01'})),...
        ScreeningData.lengte);
    HeightALL(isnan(HeightALL))=[];
    
    figure('Name','Height final')
    b1=histogram(HeightALL,'DisplayName','Data','FaceColor','#B45F66',...
            'EdgeColor','none');
    hold on
    ax = gca;
    xL=get(ax,'Ylim');
    h1=line([mHeight(1) mHeight(1)], xL, 'Color','r','LineStyle','--','LineWidth',2); % mean
    h2=line([mHeight(1)+3*sHeight(1) mHeight(1)+3*sHeight(1)], xL, 'Color','#8900BF','LineStyle','--','LineWidth',2,'MarkerEdgeColor','k'); % +ST dev
    h3=line([mHeight(1)-3*sHeight(1) mHeight(1)-3*sHeight(1)], xL, 'Color','#8900BF','LineStyle','--','LineWidth',2,'MarkerEdgeColor','k'); % +ST dev
    h4=line([mHeight(1)+4*sHeight(1) mHeight(1)+4*sHeight(1)], xL, 'Color','b','LineStyle','--','LineWidth',2); % +ST dev
    h5=line([mHeight(1)-4*sHeight(1) mHeight(1)-4*sHeight(1)], xL, 'Color','b','LineStyle','--','LineWidth',2); % +ST dev
    legend([b1 h1 h3 h5],'Real data','mean','3SD','4SD','Location','bestoutside')
    title('Height final distribution')
    xlabel('Height (cm)','FontSize',18,'FontWeight','bold') % x-axis label
    ylabel('Frequency','FontSize',18,'FontWeight','bold') % y-axis label
    ax.FontSize = 14; ax.FontWeight = 'bold';
    hold off
        
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1)= 5; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 3; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 35; %Width
    fig.Position(4) = 22; %Height
    fig.PaperOrientation = 'landscape';
        
%     print([subfold,'/',d,'/Plots/Heightfinal'],'-deps')
%     print([subfold,'/',d,'/Plots/Heightfinal','c'],'-depsc')
%     print([subfold,'/',d,'/Plots/Heightfinal'], '-dpdf','-bestfit')
    
    % 69.2725 < BMI < 18.6575
    g3 = table2array(varfun(@(x) x<(mBMI(1)-4*sBMI(1)) | x>(mBMI(1)+4*sBMI(1)),PhMsData(:,{'antwoord', 'value01', 'value02'})));
    rows_antwoord_B = idx_BMI_PhMs & g3(:,1);
    rows_value01_B = idx_BMI_PhMs & g3(:,2);
    rows_value02_B = idx_BMI_PhMs & g3(:,3);
    PhMsData.antwoord(rows_antwoord_B)=NaN;
    PhMsData.value01(rows_value01_B)=NaN;
    PhMsData.value02(rows_value02_B)=NaN;
    
    h3 = table2array(varfun(@(x) x<(mBMI(1)-4*sBMI(1)) | x>(mBMI(1)+4*sBMI(1)),PPOSdata(:,{'antwoord', 'Value01'})));
    rows_antwoord2_B = idx_BMI_PPOS & h3(:,1);
    rows_Value01_B = idx_BMI_PPOS & h3(:,2);
    PPOSdata.antwoord(rows_antwoord2_B)=NaN;
    PPOSdata.Value01(rows_Value01_B)=NaN;
    
    % Check distribution
    BMIall = vertcat(table2array(PPOSdata(idx_BMI_PPOS,{'antwoord'})),...
        table2array(PhMsData(idx_BMI_PhMs,{'antwoord'})));
    BMIall(isnan(BMIall))=[];
    BMIall(~any(BMIall,2))=[];
    
    figure('Name','BMI final')
    b1=histogram(BMIall,'DisplayName','BMI (kg/m^{2})','FaceColor','#A6E100',...
            'EdgeColor','none','FaceAlpha',0.7);
    hold on
    ax = gca;
    xL=get(ax,'Ylim');
    h1=line([mBMI(1) mBMI(1)], xL, 'Color','r','LineStyle','--','LineWidth',2); % mean
    h2=line([mBMI(1)+3*sBMI(1) mBMI(1)+3*sBMI(1)], xL, 'Color','#8900BF','LineStyle','--','LineWidth',2,'MarkerEdgeColor','k'); % +ST dev
    h3=line([mBMI(1)-3*sBMI(1) mBMI(1)-3*sBMI(1)], xL, 'Color','#8900BF','LineStyle','--','LineWidth',2,'MarkerEdgeColor','k'); % +ST dev
    h4=line([mBMI(1)+4*sBMI(1) mBMI(1)+4*sBMI(1)], xL, 'Color','b','LineStyle','--','LineWidth',2); % +ST dev
    h5=line([mBMI(1)-4*sBMI(1) mBMI(1)-4*sBMI(1)], xL, 'Color','b','LineStyle','--','LineWidth',2); % +ST dev
    legend([b1 h1 h3 h5],'Real data','mean','3SD','4SD','Location','bestoutside')
    title('BMI final distribution')
    xlabel('BMI (kg/m^{2})','FontSize',18,'FontWeight','bold') % x-axis label
    ylabel('Frequency','FontSize',18,'FontWeight','bold') % y-axis label
    ax.FontSize = 14; ax.FontWeight = 'bold';
    hold off
        
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1)= 5; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 3; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 35; %Width
    fig.Position(4) = 22; %Height
    fig.PaperOrientation = 'landscape';
        
%     print([subfold,'/',d,'/Plots/BMIfinal'],'-deps')
%     print([subfold,'/',d,'/Plots/BMIfinal','c'],'-depsc')
%     print([subfold,'/',d,'/Plots/BMIfinal'], '-dpdf','-bestfit')
    
    % Abdominal Circumference > 196.9621 cm & < 58.4506
    g4 = table2array(varfun(@(x) x<(mWC(1)-4*sWC(1)) | x>(mWC(1)+4*sWC(1)),PhMsData(:,{'value01'})));
    rows_value01_wc = idx_buik_PhMs & g4;
    PhMsData.value01(rows_value01_wc)=NaN;
    
    buikomv_80 = ScreeningData.buikomv<(mWC(1)-4*sWC(1));
    ScreeningData.buikomv(buikomv_80)= NaN;
    buikomv_200 = ScreeningData.buikomv>(mWC(1)+4*sWC(1));
    ScreeningData.buikomv(buikomv_200)= NaN;
    
    % Check ditribution
    WCall = table2array(PhMsData(idx_buik_PhMs,{'value01'}));
    WCall(isnan(WCall))=[];
    WCall(~any(WCall,2))=[];
    
    figure('Name','Height final')
    b1=histogram(WCall,'DisplayName','Data','FaceColor','#CF0047',...
            'EdgeColor','none');
    hold on
    ax = gca;
    xL=get(ax,'Ylim');
    h1=line([mWC(1) mWC(1)], xL, 'Color','r','LineStyle','--','LineWidth',2); % mean
    h2=line([mWC(1)+3*sWC(1) mWC(1)+3*sWC(1)], xL, 'Color','#8900BF','LineStyle','--','LineWidth',2,'MarkerEdgeColor','k'); % +ST dev
    h3=line([mWC(1)-3*sWC(1) mWC(1)-3*sWC(1)], xL, 'Color','#8900BF','LineStyle','--','LineWidth',2,'MarkerEdgeColor','k'); % +ST dev
    h4=line([mWC(1)+4*sWC(1) mWC(1)+4*sWC(1)], xL, 'Color','b','LineStyle','--','LineWidth',2); % +ST dev
    h5=line([mWC(1)-4*sWC(1) mWC(1)-4*sWC(1)], xL, 'Color','b','LineStyle','--','LineWidth',2); % +ST dev
    legend([b1 h1 h3 h5],'Real data','mean','3SD','4SD','Location','bestoutside')
    title('Waist circumference final distribution')
    xlabel('Waist circumference (cm)','FontSize',18,'FontWeight','bold') % x-axis label
    ylabel('Frequency','FontSize',18,'FontWeight','bold') % y-axis label
    ax.FontSize = 14; ax.FontWeight = 'bold';
    hold off
        
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1)= 5; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 3; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 35; %Width
    fig.Position(4) = 22; %Height
    fig.PaperOrientation = 'landscape';
        
%     print([subfold,'/',d,'/Plots/WCfinal'],'-deps')
%     print([subfold,'/',d,'/Plots/WCfinal','c'],'-depsc')
%     print([subfold,'/',d,'/Plots/WCfinal'], '-dpdf','-bestfit')
    
    close all
    
    clearvars g1 h1 rows_antwoord_w g2 g3 h2 h3 rows_value02_h ...
        rows_antwoord2_B rows_antwoord2_h rows_Value01_w rows_value02_B...
        rows_antwoord2_w rows_antwoord_B rows_antwoord_h rows_antwoord_w ...
        rows_value01_B rows_Value01_B rows_value01_h rows_Value01_h ...
        rows_value01_w rows_value02_w g4 rows_value01_wc buikomv_80 ...
        gewicht_220 gewicht_70 hooggew_70 hooggew_220 fig h1 h2 h3 h4 h5 ...
        xL buikomv_200;

%% Allocate values

for i=1:length(patients_list)
    patient_ID = patients_list(i);
    
    id_patient_screen = table2array(LabPackTimesPatient(:,{'PatientCode'})) == patient_ID;
    screen_day = table2array(LabPackTimesPatient(id_patient_screen,{'ScreenDay'}));
    max_screen_day = screen_day - 15;
    min_screen_day = screen_day + 15;

    
    PreOpMean.PatientCode(i) = patient_ID;
    PreOpMin.PatientCode(i) = patient_ID;
    PreOpMax.PatientCode(i) = patient_ID;
    PreOpLast.PatientCode(i) = (patient_ID);
    PreOpMedian.PatientCode(i) = patient_ID;
    PreOpNM.PatientCode(i) = patient_ID; % Excatly during screening time
    PreOpNMPrior.PatientCode(i) = patient_ID; % Before screening
    PreOpNMTotal.PatientCode(i) = patient_ID; % Up to one year before
    
    % For estimating output variable
    weight_max.PatientCode(i) = patient_ID;
    bmi_max.PatientCode(i) = patient_ID;
    height_mean.PatientCode(i) = patient_ID;
    
    idx_patient_PPOS = [PPOSdata.PatientCode] == patient_ID;
    idx_patient_PhMs = [PhMsData.PatientCode] == patient_ID;
    idx_patient_Screening = [ScreeningData.PatientCode] == patient_ID; 
        
    %% Ectract weight values
    % _Weight_
    idW_PPOS = idx_patient_PPOS & idx_weight_PPOS;
    weight_PPOS = PPOSdata(idW_PPOS, {'antwoord', 'Value01','days_from_ok'});
        % Delete values before screening time window
        weight_PPOS(~arrayfun(@(x) x >= -365, weight_PPOS.days_from_ok),:)=[];
    
    n = size(weight_PPOS,1);
        weight_list1 = {};
        for k = 1:n
            tmp = table2array(weight_PPOS(k,1:end-1));
            tmp(isnan(tmp)) = [];
            weight_list1{k} = tmp;
        end
        
        weight1 = horzcat(weight_list1',table2cell(weight_PPOS(:,end)));
        weight1(any(cellfun(@isempty,weight1),2),:) = [];
        clearvars tmp n;
    
    idW_PhMs = idx_patient_PhMs & idx_weight_PhMs;
    weight_PhMs = PhMsData(idW_PhMs,{'antwoord', 'value01', 'value02',...
        'days_from_ok'});
    weight_PhMs(~any(arrayfun(@(x) x >= -365, weight_PhMs.days_from_ok),2),:)=[];
    
    n = size(weight_PhMs,1);
        weight_list2 = {};
        for k = 1:n
            tmp = table2array(weight_PhMs(k,1:end-1));
            tmp(isnan(tmp)) = [];
            weight_list2{k} = tmp;
        end
        
        weight2 = horzcat(weight_list2',table2cell(weight_PhMs(:,end)));
        weight2(any(cellfun(@isempty,weight2),2),:) = [];
        clearvars tmp n;
    
    weightG_Screen = ScreeningData(idx_patient_Screening,{'gewicht',...
        'days_from_ok'});
    weightG_Screen(~any(arrayfun(@(x) x >= -365, weightG_Screen.days_from_ok),2),:)=[];
    
    weightH_Screen = ScreeningData(idx_patient_Screening,{'hooggew',...
        'days_from_ok'});
    weightH_Screen(~any(arrayfun(@(x) x >= -365, weightH_Screen.days_from_ok),2),:)=[];
    
    weight_list3 = num2cell(vertcat(table2array(weightG_Screen), table2array(weightH_Screen)));
    
    weight_final = sortrows(vertcat(weight1, weight2, weight_list3),2);
    
    %% Extract BMI values
    % _BMI_
    idB_PPOS = idx_patient_PPOS & idx_BMI_PPOS;
    BMI_PPOS = PPOSdata(idB_PPOS, {'antwoord', 'Value01','days_from_ok'});
    BMI_PPOS(~any(arrayfun(@(x) x >= -365, BMI_PPOS.days_from_ok),2),:)=[];
    
    n = size(BMI_PPOS,1);
        BMI_list1 = {};
        for k = 1:n
            tmp = table2array(BMI_PPOS(k,1:end-1));
            tmp(isnan(tmp)) = [];
            BMI_list1{k} = tmp;
        end
        
        BMI1 = horzcat(BMI_list1',table2cell(BMI_PPOS(:,end)));
        BMI1(any(cellfun(@isempty,BMI1),2),:) = [];
        clearvars tmp n;
    
    idB_PhMs = idx_patient_PhMs & idx_BMI_PhMs;
    BMI_PhMs = PhMsData(idB_PhMs, {'antwoord', 'value01', 'value02',...
        'days_from_ok'});
    BMI_PhMs(~any(arrayfun(@(x) x >= -365, BMI_PhMs.days_from_ok),2),:)=[];
    
    n = size(BMI_PhMs,1);
        BMI_list2 = {};
        for k = 1:n
            tmp = table2array(BMI_PhMs(k,1:end-1));
            tmp(isnan(tmp)) = [];
            BMI_list2{k} = tmp;
        end
        
        BMI2 = horzcat(BMI_list2',table2cell(BMI_PhMs(:,end)));
        BMI2(any(cellfun(@isempty,BMI2),2),:) = [];
        clearvars tmp n;
    
    BMI_Screen = ScreeningData(idx_patient_Screening,{'bmi','days_from_ok'});
    BMI_Screen(~any(arrayfun(@(x) x >= -365, BMI_Screen.days_from_ok),2),:)=[];
    
    BMI_final = vertcat(BMI1, BMI2, table2cell(BMI_Screen));
    BMI_final = sortrows(BMI_final,2);
    
    clearvars BMI_Screen BMI1 BMI2 BMI_list2 BMI_PhMs idB_PhMs BMI_PPOS ...
        idB_PPOS
    
    %% Extract Height values
    % _Height_
    idH_PPOS = idx_patient_PPOS & idx_height_PPOS;
    height_PPOS = PPOSdata(idH_PPOS, {'antwoord', 'Value01','days_from_ok'});
    height_PPOS(~any(arrayfun(@(x) x >= -365, height_PPOS.days_from_ok),2),:)=[];
    
    n = size(height_PPOS,1);
        height_list1 = {};
        for k = 1:n
            tmp = table2array(height_PPOS(k,1:end-1));
            tmp(isnan(tmp)) = [];
            height_list1{k} = tmp;
        end
        
        height1 = horzcat(height_list1',table2cell(height_PPOS(:,end)));
        height1(any(cellfun(@isempty,height1),2),:) = [];
        clearvars tmp n;
        
    idH_PhMs = idx_patient_PhMs & idx_height_PhMs;
    height_PhMs = PhMsData(idH_PhMs, {'antwoord', 'value01', 'value02',...
        'days_from_ok'});
    height_PhMs(~any(arrayfun(@(x) x >= -365, height_PhMs .days_from_ok),2),:)=[];
    
    n = size(height_PhMs,1);
        height_list2 = {};
        for k = 1:n
            tmp = table2array(height_PhMs(k,1:end-1));
            tmp(isnan(tmp)) = [];
            height_list2{k} = tmp;
        end
       
        height2 = horzcat(height_list2',table2cell(height_PhMs(:,end)));
        height2(any(cellfun(@isempty,height2),2),:) = [];
        clearvars tmp n;    
    
    height_Screen = ScreeningData(idx_patient_Screening,{'lengte',...
        'days_from_ok'});
    height_Screen(~any(arrayfun(@(x) x >= -365, height_Screen .days_from_ok),2),:)=[];
    
    height_final = vertcat(height1, height2, table2cell(height_Screen));
    height_final = sortrows(height_final,2);
    
    %% Extract waist circumference values
    % _Waist_
    idWC_PhMs = idx_patient_PhMs & idx_buik_PhMs;
    waist_PhMs = PhMsData(idWC_PhMs,{'antwoord', 'value01', 'value02',...
        'days_from_ok'});
    waist_PhMs(ismissing((sum(ismissing(waist_PhMs,NaN),2)),3),:)=[];
    waist_PhMs(~any(arrayfun(@(x) x >= -365, waist_PhMs.days_from_ok),2),:)=[];
    
    n = size(waist_PhMs,1);
        waist_list1 = {};
        for k = 1:n
            tmp = table2array(waist_PhMs(k,1:end-1));
            tmp(isnan(tmp)) = [];
            waist_list1{k} = tmp;
        end
            waist_list1 = waist_list1(~cellfun('isempty',waist_list1));
            
            if isempty(waist_list1) == 0
               waist1 = horzcat(waist_list1',table2cell(waist_PhMs(:,end)));
            else
               waist1 = {};
            end
            
            waist1(any(cellfun(@isempty,waist1),2),:) = [];
        
        clearvars tmp n;
        
    waist_Screen = ScreeningData(idx_patient_Screening,{'buikomv',...
        'days_from_ok'});
    waist_Screen(~any(arrayfun(@(x) x >= -365, waist_Screen.days_from_ok),2),:)=[];
    
    if isempty(waist_list1) == 1 && isempty(waist_Screen) == 0
       waist_final = table2cell(waist_Screen);
    end
    if isempty(waist_list1) == 0 && isempty(waist_Screen) == 1
       waist_final = waist1;
    end
    if isempty(waist_list1) == 0 && isempty(waist_Screen) == 0
       waist_final = vertcat(waist1, table2cell(waist_Screen));
    end
    
    waist_final = sortrows(waist_final,2);
    
    clearvars waist_Screen waist1 waist_list1 waist_PhMs idWC_PhMs
    
    % Extract sex and age
    idx_patient_GeneralData = [GeneralData.PatientNr] == patient_ID;
    age = table2array(GeneralData(idx_patient_GeneralData,'Age'));
    
    %% Build tables
    
    for j=1:length(var_names)
        var_name = var_names{j};
    
        if (sum(idx_patient_PPOS)+sum(idx_patient_PhMs)+sum(idx_patient_Screening)) > 0
            if strcmp(var_name,'Weight') == 1
                if isempty(weight_final) == 0
                        PreOpNMPrior.(var_name)(i) = sum(~any(cellfun(@(x) x >= -365 & x<=max_screen_day, weight_final(:,2)),2));
                        PreOpNM.(var_name)(i)= sum(cellfun(@(x) x >= max_screen_day & x <= min_screen_day , weight_final(:,2))); % Before screening
                        PreOpNMTotal.(var_name)(i) = length(weight_final);
                        
                        a = cell2mat(weight_final(:,1));
                        % For estimating output variable
                        %weight_max.Max(i) = max(a);
                        weight_max.Max(i) = min(a);
                        PreOpMin.(var_name)(i) = min(a);
                        PreOpMax.(var_name)(i) = max(a);
    
                        % Standarize values
                        s = ((mWeight(1,1)+2*sWeight(1,1)) - (mWeight(1,1)-2*sWeight(1,1)))/4; % sigma = (RefUp - RefLow)/4
                        m = mWeight(1,1);        % mu = RefLow + 2s
                        for pp=1:size(a,1)
                            a(pp,1) = (a(pp,1)-m)/s;
                        end
                        clearvars s pp
                        
                        PreOpMean.(var_name)(i) = mean(a);
                        PreOpMedian.(var_name)(i) = median(a);
                        PreOpLast.(var_name)(i) = a(end,1);
                    else
                        PreOpLast.(var_name)(i)  = NaN;
                        PreOpNM.(var_name)(i) = 0;
                        PreOpMean.(var_name)(i) = NaN;
                        PreOpMin.(var_name)(i) = NaN;
                        PreOpMax.(var_name)(i) = NaN;
                        PreOpMedian.(var_name)(i) = NaN;
                        PreOpNMPrior.(var_name)(i) = 0;
                        PreOpNMTotal.(var_name)(i) = 0;
                        % For estimating output variable
                        weight_max.Max(i) = NaN;
                end
            end
            
            if strcmp(var_name,'BMI') == 1
                if isempty(BMI_final) == 0
                    PreOpNM.(var_name)(i) = sum(cellfun(@(x) x >= max_screen_day & x <= min_screen_day , BMI_final(:,2)));
                    PreOpNMPrior.(var_name)(i) = sum(~any(cellfun(@(x) x >= -365 & x <= max_screen_day, BMI_final(:,2)),2));
                    PreOpNMTotal.(var_name)(i) = length(BMI_final);
                    
                    b=cell2mat(BMI_final(:,1));
                    bmi_max.Max(i) = max(b);
                    
                    PreOpMin.(var_name)(i) = min(b);
                    PreOpMax.(var_name)(i) = bmi_max.Max(i);
                    
                    % Standarize values
                    s = ((mBMI(1,1)+2*sBMI(1,1)) - (mBMI(1,1)-2*sBMI(1,1)))/4; % sigma = (RefUp - RefLow)/4
                    m = mBMI(1,1);        % mu = RefLow + 4s
                    for pp=1:size(b,1)
                        b(pp,1) = (b(pp,1)-m)/s;
                    end
                    clearvars s m pp
                    
                    PreOpMean.(var_name)(i) = mean(b);
                    PreOpMedian.(var_name)(i) = median(b);
                    PreOpLast.(var_name)(i) = b(end,1);
                else
                    PreOpLast.(var_name)(i) = NaN;
                    PreOpMean.(var_name)(i) = NaN;
                    PreOpMin.(var_name)(i) = NaN;
                    PreOpMax.(var_name)(i) = NaN;
                    PreOpNM.(var_name)(i) = 0;
                    PreOpMedian.(var_name)(i) = NaN;
                    PreOpNMPrior.(var_name)(i) = 0;
                    PreOpNMTotal.(var_name)(i) = 0;
                    % For estimating output variable
                    bmi_max.Max(i) = NaN;
                end
            end
            
            if strcmp(var_name,'Height') == 1
                if isempty(height_final) == 0
                    PreOpNM.(var_name)(i) = sum(cellfun(@(x) x >= max_screen_day & x <= min_screen_day , height_final(:,2)));
                    PreOpNMPrior.(var_name)(i) = sum(~any(cellfun(@(x) x >= -365 & x <= max_screen_day, height_final(:,2)),2));
                    PreOpNMTotal.(var_name)(i) = length(height_final);
                    
                    try
                        c=cell2mat(height_final(:,1));
                        c=c(c<300);
                    catch
                        c=[height_final{:,1}];
                        c=c(c<300);
                        c=c(c>100);
                    end
                    height_mean.Mean(i) = mean(c,'omitnan');
                    PreOpMin.(var_name)(i) = min(c);
                    PreOpMax.(var_name)(i) = max(c);
                    
                    % Standarize values
                    s = ((mHeight(1,1)+2*sHeight(1,1)) - (mHeight(1,1)-2*sHeight(1,1)))/4; % sigma = (RefUp - RefLow)/4
                    m = mHeight(1,1);        % mu = RefLow + 4s
                    for pp=1:size(c,1)
                        c(pp,1) = (c(pp,1)-m)/s;
                    end
                    clearvars s m pp
                    
                    PreOpMean.(var_name)(i) = mean(c);
                    PreOpMedian.(var_name)(i) = median(c);
                    PreOpLast.(var_name)(i) = c(end,1);
                else
                    PreOpLast.(var_name)(i) = NaN;
                    PreOpMean.(var_name)(i) = NaN;
                    PreOpMin.(var_name)(i) = NaN;
                    PreOpMax.(var_name)(i) = NaN;
                    PreOpNM.(var_name)(i) = 0;
                    PreOpMedian.(var_name)(i) = NaN;
                    PreOpNMPrior.(var_name)(i) = 0;
                    PreOpNMTotal.(var_name)(i) = 0;
                    height_mean.Mean(i) = NaN;
                end
            end
            
            if strcmp(var_name,'WaistCircum') == 1
                dWC=cell2mat(waist_final);
                dWC(isnan(dWC),:)=[];
                
                if isempty(dWC) == 0
                    PreOpMin.(var_name)(i) = min(dWC(:,1));
                    PreOpMax.(var_name)(i) = max(dWC(:,1));
                end
                
                % Standarize values
                s = ((mWC(1,1)+2*sWC(1,1)) - (mWC(1,1)-2*sWC(1,1)))/4; % sigma = (RefUp - RefLow)/4
                m = mWC(1,1);        % mu = RefLow + 4s
                for pp=1:size(dWC,1)
                    dWC(pp,1) = (dWC(pp,1)-m)/s;
                end
                clearvars s m pp
                
                if isempty(dWC) == 0
                    PreOpLast.WaistCircum(i) = dWC(end,1);
                    PreOpNM.(var_name)(i) = sum(arrayfun(@(x) x >= max_screen_day & x <= min_screen_day , dWC(:,2)));
                    PreOpNMPrior.(var_name)(i) = sum(~any(arrayfun(@(x) x >= -365 & x <= max_screen_day, dWC(:,2)),2));
                    PreOpNMTotal.(var_name)(i) = size(dWC,1);
                    
                    PreOpMean.(var_name)(i) = mean(dWC(:,1));
                    PreOpMedian.(var_name)(i) = median(dWC(:,1));
                else
                    PreOpLast.(var_name)(i) = NaN;
                    PreOpMean.(var_name)(i) = NaN;
                    PreOpMin.(var_name)(i) = NaN;
                    PreOpMax.(var_name)(i) = NaN;
                    PreOpNM.(var_name)(i) = 0;
                    PreOpMedian.(var_name)(i) = NaN;
                    PreOpNMPrior.(var_name)(i) = 0;
                    PreOpNMTotal.(var_name)(i) = 0;
                end
            end
            
            clearvars a b c dWC;
            
            if strcmp(var_name,'BRI') == 1
                if isempty(waist_final) == 0
                    if isempty(height_final) == 0
                        % Formula proposed by Thomas et al. (2013)
                        % BRI is a number ranging from 1 (long lean body shape)...
                        % to over 16 (round more circular body shape - OBESE). 
                        waist_BRI = cell2mat(waist_final);
                        try
                            height_BRI = cell2mat(height_final);
                        catch
                            height_BRI=[height_final{:,1}];
                            height_BRI=height_BRI(height_BRI<300);
                            height_BRI=height_BRI';
                        end
                        % Assuming height has low variations
                        aBRI=(mean(height_BRI(:,1)))/2;

                        BRI_final = NaN(size(waist_BRI));
                        for x=1:size(waist_BRI,1)
                            day_waist = waist_BRI(x,2);
                            bBRI = (waist_BRI(x,1))/(2*pi);
                            cBRI = realsqrt((aBRI^2)-(bBRI^2)); % Eccentricity
                            BRI_final(x,1) = 364.2-(365.5*(cBRI/aBRI));
                            BRI_final(x,2) = day_waist;
                        end

                        clearvars bBRI cBRI aBRI x day_waist;

                        BRI_final(isnan(BRI_final(:,1)),:)=[];
                        BRI_final=sortrows(BRI_final,2);
                        
                        if isempty(BRI_final) == 0
                            PreOpLast.(var_name)(i) = BRI_final(end,1);
                            PreOpMean.(var_name)(i) = mean(BRI_final(:,1),'omitnan');
                            PreOpMin.(var_name)(i) = min(BRI_final(:,1));
                            PreOpMax.(var_name)(i) = max(BRI_final(:,1));
                            PreOpNMTotal.(var_name)(i) = size(BRI_final,1);
                            PreOpMedian.(var_name)(i) = median(BRI_final(:,1),'omitnan');
                            PreOpNM.(var_name)(i) = NaN;
                            PreOpNMPrior.(var_name)(i) = NaN;
                            
                        else
                            PreOpLast.(var_name)(i) = NaN;
                            PreOpMean.(var_name)(i) = NaN;
                            PreOpMin.(var_name)(i) = NaN;
                            PreOpMax.(var_name)(i) = NaN;
                            PreOpNMTotal.(var_name)(i) = 0;
                            PreOpMedian.(var_name)(i) = NaN;
                            PreOpNM.(var_name)(i) = NaN;
                            PreOpNMPrior.(var_name)(i) = NaN;
                        end
                    else
                        PreOpLast.(var_name)(i) = NaN;
                        PreOpMean.(var_name)(i) = NaN;
                        PreOpMin.(var_name)(i) = NaN;
                        PreOpMax.(var_name)(i) = NaN;
                        PreOpNMTotal.(var_name)(i) = 0;
                        PreOpMedian.(var_name)(i) = NaN;
                        PreOpNM.(var_name)(i) = NaN;
                        PreOpNMPrior.(var_name)(i) = NaN;
                    end
                else
                    PreOpLast.(var_name)(i) = NaN;
                    PreOpMean.(var_name)(i) = NaN;
                    PreOpMin.(var_name)(i) = NaN;
                    PreOpMax.(var_name)(i) = NaN;
                    PreOpNMTotal.(var_name)(i) = 0;
                    PreOpMedian.(var_name)(i) = NaN;
                    PreOpNM.(var_name)(i) = NaN;
                    PreOpNMPrior.(var_name)(i) = NaN;
                end
            end
            
            clearvars BRI_final;
            
            if strcmp(var_name,'ABSI') == 1
                if isempty(waist_final) == 0
                    if isempty(height_final) == 0
                        if isempty(BMI_final) == 0
                            % Formula proposed by Krakauer et al. (2012)
                            % The higher ABSI, the higher the risk of mortality ...
                            % Assuming height has low variations
                            % Using waist_BRI and height_BRI (in meters)
                            height_ABSI = realsqrt((mean(height_BRI(:,1)))*0.01);
                            BMI_ABSI = cell2mat(BMI_final);
                            BMI_ABSI(:,1) = BMI_ABSI(:,1).^(2/3);
                             
                            ABSI_final = NaN(size(BMI_ABSI));
                            for z=1:size(BMI_ABSI,1)
                                day_BMI = BMI_ABSI(z,2);
                                BMI_power = BMI_ABSI(z,1);
                                idx_day_waist = waist_BRI(:,2) == day_BMI;
                                waist_ABSI = (waist_BRI(idx_day_waist,1))*0.01;
                                
                                if isempty(waist_ABSI) == 0
                                    ABSI_final(z,1) = mean(waist_ABSI)/(BMI_power*height_ABSI);
                                    ABSI_final(z,2) = day_BMI;
                                else
                                    ABSI_final(z,1) = NaN;
                                    ABSI_final(z,2) = day_BMI;
                                end
                            end
                    
                            clearvars waist_ABSI BMI_power idx_day_waist ...
                                day_BMI z height_ABSI BMI_ABSI;
                    
                            ABSI_final(isnan(ABSI_final(:,1)),:)=[];
                            ABSI_final=sortrows(ABSI_final,2);
                            
                            if isempty(ABSI_final) == 0
                                PreOpLast.(var_name)(i) = ABSI_final(end,1);
                                PreOpMean.(var_name)(i) = mean(ABSI_final(:,1),'omitnan');
                                PreOpMin.(var_name)(i) = min(ABSI_final(:,1));
                                PreOpMax.(var_name)(i) = max(ABSI_final(:,1));
                                PreOpNM.(var_name)(i) = size(ABSI_final,1);
                                PreOpMedian.(var_name)(i) = median(ABSI_final(:,1),'omitnan');
                                PreOpNM.(var_name)(i) = NaN;
                                PreOpNMPrior.(var_name)(i) = NaN;
                            else
                                PreOpLast.(var_name)(i) = NaN;
                                PreOpMean.(var_name)(i) = NaN;
                                PreOpMin.(var_name)(i) = NaN;
                                PreOpMax.(var_name)(i) = NaN;
                                PreOpNMTotal.(var_name)(i) = 0;
                                PreOpMedian.(var_name)(i) = NaN;
                                PreOpNM.(var_name)(i) = NaN;
                                PreOpNMPrior.(var_name)(i) = NaN;
                            end
                        else
                            PreOpLast.(var_name)(i) = NaN;
                            PreOpMean.(var_name)(i) = NaN;
                            PreOpMin.(var_name)(i) = NaN;
                            PreOpMax.(var_name)(i) = NaN;
                            PreOpNMTotal.(var_name)(i) = 0;
                            PreOpMedian.(var_name)(i) = NaN;
                            PreOpNM.(var_name)(i) = NaN;
                            PreOpNMPrior.(var_name)(i) = NaN;
                        end
                    else
                        PreOpLast.(var_name)(i) = NaN;
                        PreOpMean.(var_name)(i) = NaN;
                        PreOpMin.(var_name)(i) = NaN;
                        PreOpMax.(var_name)(i) = NaN;
                        PreOpNMTotal.(var_name)(i) = 0;
                        PreOpMedian.(var_name)(i) = NaN;
                        PreOpNM.(var_name)(i) = NaN;
                        PreOpNMPrior.(var_name)(i) = NaN;
                    end
                else
                    PreOpLast.(var_name)(i) = NaN;
                    PreOpMean.(var_name)(i) = NaN;
                    PreOpMin.(var_name)(i) = NaN;
                    PreOpMax.(var_name)(i) = NaN;
                    PreOpNMTotal.(var_name)(i) = 0;
                    PreOpMedian.(var_name)(i) = NaN;
                    PreOpNM.(var_name)(i) = NaN;
                    PreOpNMPrior.(var_name)(i) = NaN;
                end
            end
            
            clearvars ABSI_final;
                                    
            if strcmp(var_name,'TBFM') == 1
                if isempty(BMI_final) == 0
                    if isempty(weight_final) == 0
                       % Formula proposed by DEURENBERG et al. (2006)
                       % Adults: (BF?%)?=?(1.2 × BMI)?+?(0.23 × age)???5.4
                       BMI_TBFM = cell2mat(BMI_final);
                       weight_TBFM = cell2mat(weight_final);
                                                    
                       TBFM_final = NaN(size(BMI_TBFM));
                        for w=1:size(BMI_TBFM,1)
                            day_BMI = BMI_TBFM(w,2);
                            BMI_BF = BMI_TBFM(w,1);
                            idx_day_weight = weight_TBFM(:,2) == day_BMI;
                            weight_BF = mean(weight_TBFM(idx_day_weight,1),'omitnan');

                            if isempty(BMI_BF) == 0
                                if isempty(weight_BF) == 0
                                    TBFM_final(w,1) = (((1.2*BMI_BF)+(0.23*age)-5.4)*0.01)*weight_BF;
                                    TBFM_final(w,2) = day_BMI;
                                else
                                    TBFM_final(w,1) = NaN;
                                    TBFM_final(w,2) = day_BMI;
                                end                             
                            else
                                TBFM_final(w,1) = NaN;
                                TBFM_final(w,2) = day_BMI;
                            end
                        end
                    
                        clearvars BMI_BF idx_day_weight weight_BF day_BMI w;
                    
                        TBFM_final(isnan(TBFM_final(:,1)),:)=[];
                        TBFM_final=sortrows(TBFM_final,2);
                            
                        if isempty(TBFM_final) == 0
                            PreOpLast.(var_name)(i) =  TBFM_final(end,1);
                            PreOpMean.(var_name)(i) = mean(TBFM_final(:,1),'omitnan');
                            PreOpMin.(var_name)(i) = min(TBFM_final(:,1));
                            PreOpMax.(var_name)(i) = max(TBFM_final(:,1));
                            PreOpNMTotal.(var_name)(i) = size(TBFM_final,1);
                            PreOpMedian.(var_name)(i) = median(TBFM_final(:,1),'omitnan');
                            PreOpNM.(var_name)(i) = NaN;
                            PreOpNMPrior.(var_name)(i) = NaN;
                        else
                            PreOpLast.(var_name)(i) = NaN;
                            PreOpMean.(var_name)(i) = NaN;
                            PreOpMin.(var_name)(i) = NaN;
                            PreOpMax.(var_name)(i) = NaN;
                            PreOpNMTotal.(var_name)(i) = 0;
                            PreOpMedian.(var_name)(i) = NaN;
                            PreOpNM.(var_name)(i) = NaN;
                            PreOpNMPrior.(var_name)(i) = NaN;
                        end
                    else
                        PreOpLast.(var_name)(i) = NaN;
                        PreOpMean.(var_name)(i) = NaN;
                        PreOpMin.(var_name)(i) = NaN;
                        PreOpMax.(var_name)(i) = NaN;
                        PreOpNMTotal.(var_name)(i) = 0;
                        PreOpMedian.(var_name)(i) = NaN;
                        PreOpNM.(var_name)(i) = NaN;
                        PreOpNMPrior.(var_name)(i) = NaN;
                    end
                else
                    PreOpLast.(var_name)(i) = NaN;
                    PreOpMean.(var_name)(i) = NaN;
                    PreOpMin.(var_name)(i) = NaN;
                    PreOpMax.(var_name)(i) = NaN;
                    PreOpNMTotal.(var_name)(i) = 0;
                    PreOpMedian.(var_name)(i) = NaN;
                    PreOpNM.(var_name)(i) = NaN;
                    PreOpNMPrior.(var_name)(i) = NaN;
                end
            end
            
            clearvars TBFM_final;
            
        else
            PreOpLast.(var_name)(i) = NaN;
            PreOpMean.(var_name)(i) = NaN;
            PreOpMin.(var_name)(i) = NaN;
            PreOpMax.(var_name)(i) = NaN;
            PreOpNMTotal.(var_name)(i) = 0;
            PreOpMedian.(var_name)(i) = NaN;
            PreOpNM.(var_name)(i) = NaN;
            PreOpNMPrior.(var_name)(i) = NaN;
        end
    end   
end

clearvars var_name idx_patient_PPOS idx_patient_PhMs idx_patient_Screening ...
    patient_ID idx_patient_Screening height_final BMI_final waist_final ...
    weight_final idx_patient_PPOS idx_BMI_PPOS height_Screen height_list2 ...
    height2 idH_PhMs height_PhMs height_PPOS idH_PPOS height_list1 i k j ...
    height1 BMI_list1 waist1 waist_list1 waist_PhMs waist_Screen weight1 ...
    weight2 weight_list1 weight_list2 weight_list3 weight_PhMs weight_PPOS ...
    weightG_Screen weightH_Screen var_names col_names prealloc_num ...
    PhMsData waist_BRI height_BRI age BMI_TBFM weight_TBFM GeneralData ...
    ScreeningData PPOSdata;

%% Hypthesis testing

for p = 2:(length(PreOpMin.Properties.VariableNames))
    var_name = PreOpMin.Properties.VariableNames{p};
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((PreOpMin.(var_name)),...
        (PreOpMax.(var_name)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)';
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((PreOpMin.(var_name)),...
        (PreOpMax.(var_name)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)';
    else
        veredict2 = 'No diff sig (var)';
    end
    
    stat.(var_name).tstudent.fieldName = (var_name);
    stat.(var_name).tstudent.h_ttest = httest2;
    stat.(var_name).tstudent.pval_ttest = pttest2;
    stat.(var_name).tstudent.ci_ttest = cittest2;
    stat.(var_name).tstudent.tstat = statsttest2.tstat;
    stat.(var_name).tstudent.df = statsttest2.df;
    stat.(var_name).tstudent.sd = statsttest2.sd;
    stat.(var_name).tstudent.tstudent = veredict1;
    
    stat.(var_name).Ftest.fieldName = var_name;
    stat.(var_name).Ftest.fstat = stats_var.fstat;
    stat.(var_name).Ftest.df1 = stats_var.df1;
    stat.(var_name).Ftest.df2 = stats_var.df2;
    stat.(var_name).Ftest.h_ftest = h_var;
    stat.(var_name).Ftest.pval_ftest = p_var;
    stat.(var_name).Ftest.ci_ftest = ci_var;
    stat.(var_name).Ftest.ftest = veredict2;
end