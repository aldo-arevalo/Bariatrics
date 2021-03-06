function [labels] = getVarsLabels(variable_names, Options)
% Funtion entitled to rename the variable names of CZE dataset
% The hospital has assigned an unique code to each variable. However, for
% publication purposes these features names should not be used and it is
% recommended to use generic names.

% [labels] = M2_arrange_labels(variable_names)
%     Input:
%       variable_names  ... M cell vector,
%                           M is the number of features
%       Options ... Structure which contains the options selected to build
%                   the models
%
%     Output:
%       labels  ... M cell vector,
%                   M is the number of features with generic names
%
% Dependencies: getDataset > getPartitions
%
% Author: Aldo Arevalo
% Date: 10/11/2020

% Create string vectors
   labels1 = {'p[Albumin]','p[ALP]','p[ALAT]','p[ASAT]','p[Bili]',...
       'p[Cholesterol]','p[Ch_HDL]','p[C_pept]','p[CRP]','Erythrocytes',...
       'p[Ferritine]','p[Folate]','p[(PO_4)^3-]','p[Glucose]','p[GGT]',...
       'p[HDL]','Hematocrit','p[Hb]','p[HbA1c]','p[Fe]','p[Insulin]',...
       'p[K^+]','p[Creatinine]','p[LD]','p[Leukos]','p[Mg^2]','MCH',...
       'MCV','p[Na^+]','p[PTH]','p[Triglycerides]','Thrombocytes',...
       'p[TSH]','p[Urea]','p[VitB_6]','p[VitB_{12}]','p[VitD]',...
       'p[VitB_1]','p[VitA]','p[Zn]','p[Ca^2+]','GFR','Systolic','Age',...
       'Gender','Smoker','Cigarettes','Diabetes','Hypertension',...
       'Dyslipidemia','BMI','TBFM','Diastolic','vID120848','vID120881',...
       'vID120883','vID120836','hypmedn','insuleh','gewrklacht',....
       'Height','WeightMax','BMImax','WCmax','BRImax','ABSImax',...
       'TBFMmax','WeightMin','BMImin','WCmin','BRImin','ABSImin',...
       'TBFMmin','Comorb_Score','EtOH_Score','QoL_score','BARO_score',...
       'Rand36_score','vID120838','vID120843','vID120847','vID121223',...
       'oraaln','insuln','dyslipmedn','p[Cortisol]','p[LDL]'};
   
   labels2 = {'BAL002','BAL014','BAL015','BAS002','BBI001','BCH004',...
       'BCH014','BCP002','BCR002','BER002','BFE006','BFO000','BFO002',...
       'BGL003','BGT001','BHD000','BHE000','BHE001','BHE030','BIJ001',...
       'BIN000','BKA000','BKR000','BLD004','BLE001','BMA002','BMC000',...
       'BMC002','BNA001','BPA016','BTR003','BTR006','BTS003','BUR002',...
       'BVI004','BVI028','BVI030','BVI031','BVI032','BZI003','BCA006b',...
       'BKR000b','Systolic','Age','Sex','roken','rokenpy','diabet',...
       'hypert','dyslip','BMI','TBFM','Diastolic','Hist_Hyper',...
       'Hist_Other','Hist_Dr','Hist_Surg','Hypertension_Med',...
       'Insulin_Units','Arthros','Height','WeightMax','BMImax','WCmax',...
       'BRImax','ABSImax','TBFMmax','WeightMin','BMImin','WCmin',...
       'BRImin','ABSImin','TBFMmin','Comorb_Score','EtOH_Score',...
       'QoL_score','BARO_score','Rand36_score','Hist_Second','Hist_Chest',...
       'Hist_Heart','Hist_Dent','Antidiab_Med','Insulin_Med',...
       'Dyslip_Med','BCO017','BLD002'};	
   
   if (length(labels1) == length(labels2)) == true
       % Verify if it contains lables PatientCode and any metric
       idx_patient = strcmp(variable_names, 'PatientCode');
       variable_names(idx_patient) = [];
       idx_patient = strcmp(variable_names, 'dyslip');
       variable_names(idx_patient) = [];
       idx_patient = strcmp(variable_names, 'hypert');
       variable_names(idx_patient) = [];
       idx_patient = strcmp(variable_names, 'Age_');
       variable_names(idx_patient) = [];

        if Options.TWL == true
          idx_metric = strcmp(variable_names, 'TWL');
        end
        if Options.EBWL == true
          idx_metric = strcmp(variable_names, 'EWL');
        end
        if Options.EBMIL == true
          idx_metric = strcmp(variable_names, 'EBMIL');
        end
        if Options.WL == true
          idx_metric = strcmp(variable_names, 'WL');
        end

       variable_names(idx_metric) = [];

       % Compare the input againts all_labels
       labels = {};
       for p=1:length(variable_names)
           var_name = variable_names{p};
           idx_var = strcmp(labels2,var_name);
           labels{1,p} = labels1{idx_var};
       end
   
   end
   