%% Building the table containing all analysis per patient
% Obj: Make unique list U with BEP codes
% Minimal frequency >=100 Instances. Treatment done in SQL

load([subfold,'/',d,'/LabZt.mat'],'LabResults', 'GeneralData',...
    'Calcium','ASAT','Albumin','ALAT','BKR00b', 'SD')

%% Initialize Structure
% Extract bepcodes in a cell array
var_names=unique(LabResults.bepcode);

data_struct.Result = {};
data_struct.days_from_ok = NaN;

ColumSruct.PatientCode = NaN;
for i = 1:length(var_names)
    ColumSruct.(var_names{i}) = data_struct;
end

%% Fill in the structure for each patient

list_PatientID = unique(LabResults.PatientCode);

% Pre-allocate structs for all patients
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
        
        vars = table2array(s(idx_var,{'StandardizedResult'}));
        %vars = table2array(s(idx_var,{'uitslagNUM'}));
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

%% Transform lab analyses: sCr, Free [sCa2+], estimation of GFR, INR, ASAT/ALAT

% CKD-EPI constants (Levy et al, 2009)
   % Assuming all patients are caucasians
   kappa_female = 0.7;
   kappa_male = 0.9;
   alpha_female = -0.329;
   alpha_male = -0.411;
   betha_female = 1.018;
   betha_male = 1;
   
   % Conversion factor to obtain [sCr] in mg/dL from uM
   conv_factor = (1/1e6)*(113.2)*(1000)*(1/10);
   
for i=1:length(list_PatientID)
    patient_ID = LabResultsPerPatient(i).PatientCode;
    
    % Calculation of Free serum Calcium
    idx_patient_Alb_cel = [Albumin.PatientCode] == patient_ID;
    idx_patient_Cal_cel = [Calcium.PatientCode] == patient_ID;
    
    Alb_meas = table2array(Albumin(idx_patient_Alb_cel,{'uitslagNUM','days_from_ok'}));
    Alb_reflow = table2array(Albumin(idx_patient_Alb_cel,{'RefLow'}));
    Alb_refup = table2array(Albumin(idx_patient_Alb_cel,{'RefUp'}));
    
    sCa = table2array(Calcium(idx_patient_Cal_cel,{'uitslagNUM','days_from_ok'}));
    Cal_reflow = table2array(Calcium(idx_patient_Cal_cel,{'RefLow'}));
    Cal_refup = table2array(Calcium(idx_patient_Cal_cel,{'RefUp'}));
    
        if isempty(sCa) == 0
            if isempty(Alb_meas) == 0                
               free_calcium = NaN(size(sCa,1),3);
                    for j=1:size(sCa,1)
                        day_calcium = sCa(j,2);
                        total_Calcium = sCa(j,1); % mM
                        total_Calcium_low = Cal_reflow(j,1); %mM
                        total_Calcium_up = Cal_refup(j,1); % mM
                        idx_day_albumin = Alb_meas(:,2) == day_calcium;
                        
                        if sum(idx_day_albumin) == 0
                            free_calcium(j,1) = NaN;
                            free_calcium(j,2) = NaN;
                            free_calcium(j,3) = NaN;
                        else
                            total_albumin = mean(Alb_meas(idx_day_albumin,1)); % g/L
                            total_albumin_up = mean(Alb_refup(idx_day_albumin,1)); % g/L
                            total_albumin_down = mean(Alb_reflow(idx_day_albumin,1)); % g/L
                            
                            free_calcium(j,1) = total_Calcium + (0.02*(40-total_albumin));
                            free_calcium(j,2) = total_Calcium_low + (0.02*(40-total_albumin_down));
                            free_calcium(j,3) = total_Calcium_up + (0.02*(40-total_albumin_up));
                        end
                    end
                    clearvars j total_albumin total_albumin_up ...
                        total_albumin_down total_Calcium_low ...
                        total_Calcium_up
            
                % Standarize values (Z-transformation)
                zsCa = NaN(size(free_calcium,1),1);

                for mm=1:size(free_calcium,1)
                    if and(~isnan(free_calcium(mm,2)),~isnan(free_calcium(mm,3)))
                        % Both values are available
                        s = (free_calcium(mm,3) - free_calcium(mm,2))/4; % sigma = (RefUp - RefLow)/4
                        m = free_calcium(mm,2) + SD*s;            % mu = RefLow + 2s
                        zsCa(mm,1) = (free_calcium(mm,1)-m)/s;
                    elseif and(~isnan(free_calcium(mm,2)),isnan(free_calcium(mm,3)))
                        % RefUp not available, one sided/tail 
                        zsCa(mm,1) = (free_calcium(mm,1)/free_calcium(mm,2))*4-6;
                    end
                end
            end
            
            clearvars mm free_calcium
            
            LabResultsPerPatient(i).BCA006b.Result = zsCa;
            LabResultsPerPatient(i).BCA006b.days_from_ok = sCa(:,2);

        else
            LabResultsPerPatient(i).BCA006b.Result = NaN;
            LabResultsPerPatient(i).BCA006b.days_from_ok = NaN;
        end
    clearvars j day_calcium idx_day_albumin sCa zsCa Alb_reflow ...
       total_albumin Alb_meas Alb_refup Cal_reflow Cal_refup ...
       idx_patient_Alb_cel idx_patient_Cal_cel m s total_Calcium ...
       total_Calcium_low total_Calcium_up
   
   % Estimation of ASAT/ALAT ratio
    idx_patient_ASAT = [ASAT.PatientCode] == patient_ID;
    idx_patient_ALAT = [ALAT.PatientCode] == patient_ID;
    
    ASAT_m = table2array(ASAT(idx_patient_ASAT,{'uitslagNUM','days_from_ok'}));
    ASAT_reflow = table2array(ASAT(idx_patient_ASAT,{'RefLow'}));
    ASAT_refup = table2array(ASAT(idx_patient_ASAT,{'RefUp'}));
    
    ALAT_m = table2array(ALAT(idx_patient_ALAT,{'uitslagNUM','days_from_ok'}));
    ALAT_refup = table2array(ALAT(idx_patient_ALAT,{'RefUp'}));
    
        if isempty(ALAT_m) == 0
            if isempty(ASAT_m) == 0               
                ASAT_ALATr = NaN(size(ALAT_m,1),3);
                
                for p=1:size(ALAT_m,1)
                    day_ALAT = ALAT_m(p,2);
                    ALAT_meas = ALAT_m(p,1); % U/L
                    ALAT_up = ALAT_refup(p,1); % U/L
                    idx_day_ASAT = ASAT_m(:,2) == day_ALAT;
                    
                    if sum(idx_day_ASAT) == 0
                        ASAT_ALATr(p,1) = NaN;
                        ASAT_ALATr(p,2) = NaN;
                        ASAT_ALATr(p,3) = NaN;
                    else
                        ASAT_meas = mean(ASAT_m(idx_day_ASAT,1));
                        ASAT_low = mean(ASAT_reflow(idx_day_ASAT,1));
                        ASAT_up = mean(ASAT_refup(idx_day_ASAT,1));
                        
                        ASAT_ALATr(p,1) = ASAT_meas/ALAT_meas;
                        ASAT_ALATr(p,2) = ASAT_low;
                        ASAT_ALATr(p,3) = ASAT_up/ALAT_up;
                    end
                end
                clearvars p day_ALAT ASAT_meas ALAT_meas ALAT_up ...
                    ASAT_up idx_day_ASAT
                
                % Standarize values (Z-transformation)
                zASAT_ALAT = NaN(size(ASAT_ALATr,1),1);

                for mm=1:size(ASAT_ALATr,1)
                    if and(~isnan(ASAT_ALATr(mm,2)),~isnan(ASAT_ALATr(mm,3)))
                        % Both values are available
                        s = (ASAT_ALATr(mm,3) - ASAT_ALATr(mm,2))/4; % sigma = (RefUp - RefLow)/4
                        m = ASAT_ALATr(mm,2) + SD*s;            % mu = RefLow + 2s
                        zASAT_ALAT(mm,1) = (ASAT_ALATr(mm,1)-m)/s;
                    elseif and(~isnan(ASAT_ALATr(mm,2)),isnan(ASAT_ALATr(mm,3)))
                        % RefUp not available, one sided/tail 
                        zASAT_ALAT(mm,1) = (ASAT_ALATr(mm,1)/ASAT_ALATr(mm,2))*4-6;
                    end
                end
                
                LabResultsPerPatient(i).ASAT_ALAT.Result =  zASAT_ALAT;
                LabResultsPerPatient(i).ASAT_ALAT.days_from_ok = ALAT_m(:,2);
            else
                LabResultsPerPatient(i).ASAT_ALAT.Result =  NaN;
                LabResultsPerPatient(i).ASAT_ALAT.days_from_ok = NaN;
            end
        else
            LabResultsPerPatient(i).ASAT_ALAT.Result =  NaN;
            LabResultsPerPatient(i).ASAT_ALAT.days_from_ok = NaN;
        end
        
   clearvars zASAT_ALAT mm m ASAT_refup ASAT_reflow ASAT_m ASAT_low ...
             ASAT_ALATr ALAT_refup ALAT_reflow ALAT_m ALAT_low
   
   % Glomerular Filtration Rate (GFR) to correct sCr values
   % Forumla: CKD-EPI (Levy et al, 2009)
   % Sex: 0/false = female; 1/true = male
   % GFR = [mL/min/1.73 m2]
   
   idx_patient_GeneralData = [GeneralData.PatientNr] == patient_ID;
   idx_patient_BKR00b = [BKR00b.PatientCode] == patient_ID;
   age = table2array(GeneralData(idx_patient_GeneralData,'Age'));
   sex = GeneralData(idx_patient_GeneralData,'Sex');
   
   Creatinine_meas = table2array(BKR00b(idx_patient_BKR00b,{'uitslagNUM','days_from_ok'}));
   Creatinine_reflow = table2array(BKR00b(idx_patient_BKR00b,{'RefLow'}));
   Creatinine_refup = table2array(BKR00b(idx_patient_BKR00b,{'RefUp'}));
   
        if isempty(Creatinine_meas) == 0
            
            sCr_meas = Creatinine_meas(:,1);    
            GFR = NaN(size(sCr_meas,1),3);
                if sex{1,1} == 0 % 0=female
                    for z=1:size(sCr_meas,1)
                        sCr = (sCr_meas(z,1))*(conv_factor);
                        GFR(z,1)= 141*(min(sCr/kappa_female,1))^alpha_female*(max(sCr/kappa_female,1)^(-1.209))*(0.993^age)*betha_female;
                        
                        sCr_low = (Creatinine_reflow(z,1))*(conv_factor);
                        GFR(z,2) = 141*(min(sCr_low/kappa_female,1))^alpha_female*(max(sCr_low/kappa_female,1)^(-1.209))*(0.993^age)*betha_female;
                        
                        sCr_up = (Creatinine_refup(z,1))*(conv_factor);
                        GFR(z,3) = 141*(min(sCr_up/kappa_female,1))^alpha_female*(max(sCr_up/kappa_female,1)^(-1.209))*(0.993^age)*betha_female;
                    end                     
                else
                    for z=1:size(sCr_meas,1)
                        sCr = (sCr_meas(z,1))*(conv_factor);
                        GFR(z,1) = 141*(min(sCr/kappa_male,1))^alpha_male*(max(sCr/kappa_male,1)^(-1.209))*(0.993^age)*betha_male;
                        
                        sCr_low = (Creatinine_reflow(z,1))*(conv_factor);
                        GFR(z,2) = 141*(min(sCr_low/kappa_male,1))^alpha_male*(max(sCr_low/kappa_male,1)^(-1.209))*(0.993^age)*betha_male;
                        
                        sCr_up = (Creatinine_refup(z,1))*(conv_factor);
                        GFR(z,3) = 141*(min(sCr_up/kappa_male,1))^alpha_male*(max(sCr_up/kappa_male,1)^(-1.209))*(0.993^age)*betha_male;
                    end
                    
                end
            
            clearvars z
            
            zGFR = NaN(size(GFR,1),1);
            
            % Standarize values
            for mm=1:size(GFR,1)
                if and(~isnan(GFR(mm,2)),~isnan(GFR(mm,3)))
                    % Both values are available
                    s = (GFR(mm,3) - GFR(mm,2))/4; % sigma = (RefUp - RefLow)/4
                    m = GFR(mm,2) + SD*s;            % mu = RefLow + 2s
                    zGFR(mm,1) = (GFR(mm,1)-m)/s;
                elseif and(~isnan(GFR(mm,2)),isnan(GFR(mm,3)))
                    % RefUp not available, one sided/tail 
                    zGFR(mm,1) = (GFR(mm,1)/GFR(mm,2))*4-6;
                end
            end
            
            LabResultsPerPatient(i).BKR000b.Result = zGFR;
            LabResultsPerPatient(i).BKR000b.days_from_ok = Creatinine_meas(:,2);
            
        else
            LabResultsPerPatient(i).BKR000b.Result =  NaN;
            LabResultsPerPatient(i).BKR000b.days_from_ok =  NaN;
        end
        
        clearvars z age sex idx_patient_GeneralData GFR sCr sCr_empty ...
            sCr_meas;
end

clearvars alpha_female alpha_male betha_female betha_male conv_factor ...
    i kappa_female kappa_male patient_ID p idx_patient_BKR00b ...
    idx_patient_ASAT idx_patient_ALAT Creatinine_meas Creatinine_reflow ...
    Creatinine_refup mm m s sCr_low sCr_up zGFR;

%%
% _Created by Aldo Arévalo_ 