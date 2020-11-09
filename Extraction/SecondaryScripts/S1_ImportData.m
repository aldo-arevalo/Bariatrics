%% Import general data from spreadsheet
% Script for importing data from the following spreadsheet:
%    Workbook: DATO CZE 2016-03.xlsx
%    Worksheet: patient
% Import the data, extracting spreadsheet dates in Excel serial date format

ver = version;
TF = contains(ver,'R2020')||contains(ver,'R2019');
if TF == 1
    %% Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 4);

    % Specify sheet and range
    opts.Sheet = "patient";
    opts.DataRange = "A5:D3613";

    % Specify column names and types
    opts.VariableNames = ["PatientNr", "OKDate", "Age", "Sex"];
    opts.SelectedVariableNames = ["PatientNr", "OKDate", "Age", "Sex"];
    opts.VariableTypes = ["double", "datetime", "double", "categorical"];
    opts = setvaropts(opts, 2, "InputFormat", "");
    opts = setvaropts(opts, 4, "EmptyFieldRule", "auto");

    % Import the data
    GeneralData = readtable("DATO CZE 2016-03.xlsx", opts, "UseExcel", false);
    
    % Clear temporary variables
    clearvars opts;
else
    [~, ~, raw, dates] = xlsread('DATO CZE 2016-03.xlsx','patient','A5:D3613','',@convertSpreadsheetExcelDates);
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,4);
    raw = raw(:,[1,3]);
    dates = dates(:,2);

    % Create output variable
    data = reshape([raw{:}],size(raw));

    % Create table
    GeneralData = table;

    % Allocate imported array to column variable names
    GeneralData.PatientNr = data(:,1);
    GeneralData.OKDate = datenum(datetime([dates{:,1}].', 'ConvertFrom', 'Excel', 'Format', 'dd/MM/yyyy'));
    GeneralData.Age = data(:,2);

    GeneralData.Sex = cellVectors(:,1);
    
    % Clear temporary variables
    clearvars data raw dates cellVectors;
end

GeneralData.Properties.VariableDescriptions{4} = 'female/vrouw=0; male/mannetje=1';
GeneralData.Sex = ismissing(GeneralData.Sex,'M');

% Save raw data
save([subfold,'/',d,'/RawData/GeneralDataRAW'],'GeneralData');

if TF == 1
    % Remove patients before 01-03-2012
    DateString1 = '01-Mar-2012';
    formatIn = 'dd-MMM-yyyy';
    Date1 = datetime(DateString1,'InputFormat',formatIn);
    GeneralData(([GeneralData.OKDate] < Date1),:)=[];
    
    % Remove patients after 01-03-2014
    DateString2 = '01-Mar-2014';
    Date2 = datetime(DateString2,'InputFormat',formatIn);
    GeneralData(([GeneralData.OKDate] > Date2),:)=[];
    
else
    % Remove patients before 01-03-2012
    DateString1 = '01-Mar-2012';
    formatIn = 'dd-mmm-yyyy';
    Date1 = datenum(DateString1,formatIn);
    GeneralData(([GeneralData.OKDate] < Date1),:)=[];

    % Remove patients after 01-03-2014
    DateString2 = '01-Mar-2014';
    Date2 = datenum(DateString2,formatIn);
    GeneralData(([GeneralData.OKDate] > Date2),:)=[];
end

% Obtain total patients list
total_patientsID = unique(GeneralData.PatientNr);

clearvars DateString1 DateString2 formatIn Date1 Date2;
%% Import screening data from spreadsheet
% Script for importing data from the following spreadsheet:
%    Original Workbook: DATO CZE 2016-03.xlsx
%    Worksheet: screening
%    Extraction Workbook: ScreeningDATO CZE2017.xlsx (contains "days from
%    ok" variable.

% Import screening data

if TF == 1
    % Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 70);

    % Specify sheet and range
    opts.Sheet = "screening";
    opts.DataRange = "A2:BR3609";

    % Specify column names and types
    opts.VariableNames = ["PatientCode", "datbez1", "gewicht", "lengte", "bmi", "buikomv", "hooggew", "roken", "rokenpy", "alcohol", "alcoholeh", "hypert", "hypmed", "hypmedn", "diabet", "diahba1c", "diamed", "oraal", "oraaln", "insul", "insuln", "insuleh", "dyslip", "dyslipldl", "dysliphdl", "dysliptrig", "dyslipratio", "dyslipmed", "dyslipmedn", "oesofa", "oesofadiag", "oesofamed", "osas", "cpap", "gewrklacht", "comorb", "comorbsetall", "comorbcar", "comorbvas", "comorbpul", "commda", "comhep", "psych", "comtro1", "comcar04", "comcar01", "comcar08", "corafw", "comcar07", "comneu1", "comvas3", "paod", "aa", "compul1", "fibro", "compul3", "hernia", "upept", "hpylor", "cirros", "galsteen", "schiz", "depri", "bipol", "angst", "menstru", "comorbove", "overigcomorb", "laparovg", "days_from_ok"];
    opts.SelectedVariableNames = ["PatientCode", "datbez1", "gewicht", "lengte", "bmi", "buikomv", "hooggew", "roken", "rokenpy", "alcohol", "alcoholeh", "hypert", "hypmed", "hypmedn", "diabet", "diahba1c", "diamed", "oraal", "oraaln", "insul", "insuln", "insuleh", "dyslip", "dyslipldl", "dysliphdl", "dysliptrig", "dyslipratio", "dyslipmed", "dyslipmedn", "oesofa", "oesofadiag", "oesofamed", "osas", "cpap", "gewrklacht", "comorb", "comorbsetall", "comorbcar", "comorbvas", "comorbpul", "commda", "comhep", "psych", "comtro1", "comcar04", "comcar01", "comcar08", "corafw", "comcar07", "comneu1", "comvas3", "paod", "aa", "compul1", "fibro", "compul3", "hernia", "upept", "hpylor", "cirros", "galsteen", "schiz", "depri", "bipol", "angst", "menstru", "comorbove", "overigcomorb", "laparovg", "days_from_ok"];
    opts.VariableTypes = ["double", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "double", "double", "double", "double", "double", "double", "string", "double", "double"];
    opts = setvaropts(opts, 2, "InputFormat", "");
    opts = setvaropts(opts, [37, 60, 61, 68], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [37, 60, 61, 68], "EmptyFieldRule", "auto");

    % Import the data
    ScreeningData = readtable("ScreeningDATO CZE2017.xlsx", opts, "UseExcel", false);
    
    % Select columns
    ScreeningData = ScreeningData(:,["PatientCode", "datbez1", "gewicht"...
        , "lengte" , "bmi", "buikomv", "hooggew", "roken", "rokenpy"...
        , "alcohol" , "alcoholeh", "hypert", "hypmed", "hypmedn"...
        , "diabet", "diahba1c" , "diamed", "oraal", "oraaln", "insul"...
        , "insuln", "insuleh", "dyslip", "dyslipldl", "dysliphdl"...
        , "dysliptrig", "dyslipratio", "dyslipmed", "dyslipmedn"...
        , "oesofa", "oesofadiag", "oesofamed", "osas", "cpap"...
        , "gewrklacht", "comorb", "comorbcar", "comorbvas", "comorbpul"...
        , "commda", "comhep", "psych", "comtro1", "comcar04", "comcar01"...
        , "comcar08", "corafw", "comcar07", "compul1", "fibro", "compul3"...
        , "schiz", "depri", "bipol", "angst", "menstru", "comorbove"...
        , "laparovg", "days_from_ok"]);
else
    % Import data
        [~, ~, raw0_0, dates0_0] = xlsread('ScreeningDATO CZE2017.xlsx','screening','A2:AJ3609','',@convertSpreadsheetExcelDates);
            % AK comborbsetall = empty
        [~, ~, raw0_1, dates0_1] = xlsread('ScreeningDATO CZE2017.xlsx','screening','AL2:AW3609','',@convertSpreadsheetExcelDates);
            % AX comneu1 = only class 1
            % AY comvas3 = only class 1
            % AZ paod = only class 1
            % BA aa = only class 1,
        [~, ~, raw0_2, dates0_2] = xlsread('ScreeningDATO CZE2017.xlsx','screening','BB2:BD3609','',@convertSpreadsheetExcelDates);
            % BE hernia = only class 1
            % BF upept = only class 1
            % BG hpylor = only class 1
            % BH cirros = empty
            % BI galsteen = empty
        [~, ~, raw0_3, dates0_3] = xlsread('ScreeningDATO CZE2017.xlsx','screening','BJ2:BO3609','',@convertSpreadsheetExcelDates);
            % BP overigcomorb = string variables
        [~, ~, raw0_4, dates0_4] = xlsread('ScreeningDATO CZE2017.xlsx','screening','BQ2:BR3609','',@convertSpreadsheetExcelDates);
        raw = [raw0_0,raw0_1,raw0_2,raw0_3,raw0_4];
        dates = [dates0_0,dates0_1,dates0_2,dates0_3,dates0_4];
        raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
        cellVectors = raw(:,1);
        raw = raw(:,[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,...
            24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,...
            46,47,48,49,50,51,52,53,54,55,56,57,58,59]);
        dates = dates(:,2);

    % Replace non-numeric cells with NaN
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
        raw(R) = {NaN}; % Replace non-numeric cells
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
        dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN

    % Create output variable
        data = reshape([raw{:}],size(raw));

    % Create table
        ScreeningData = table;
    
    % Allocate imported array to column variable names
        ScreeningData.PatientCode = str2double(cellVectors(:,1));
        ScreeningData.Properties.VariableDescriptions{1} = 'Patient ID';

        ScreeningData.datbez1 = datenum(datetime([dates{:,1}].', 'ConvertFrom', 'Excel'));
        ScreeningData.Properties.VariableDescriptions{2} = 'Date of consultation for screening';

        ScreeningData.gewicht = data(:,1);
        ScreeningData.Properties.VariableUnits{3} = 'kg';
        ScreeningData.Properties.VariableDescriptions{3} = 'Weight before surgery';

        ScreeningData.lengte = data(:,2);
        ScreeningData.Properties.VariableUnits{4} = 'm';
        ScreeningData.Properties.VariableDescriptions{4} = 'Height';

        ScreeningData.bmi = data(:,3);
        ScreeningData.Properties.VariableUnits{5} = 'kg/m^2';
        ScreeningData.Properties.VariableDescriptions{5} = 'BMI before surgery';

        ScreeningData.buikomv = data(:,4);
        ScreeningData.Properties.VariableUnits{6} = 'cm';
        ScreeningData.Properties.VariableDescriptions{6} = 'Weist circumference';

        ScreeningData.hooggew = data(:,5);
        ScreeningData.Properties.VariableUnits{7} = 'kg';
        ScreeningData.Properties.VariableDescriptions{7} = 'Greatest weight';

        ScreeningData.roken = data(:,6);
        ScreeningData.Properties.VariableDescriptions{8} = 'Smoking patient? 0=Never, 1=Yes, 2=Stopped';

        ScreeningData.rokenpy = data(:,7);
        ScreeningData.Properties.VariableUnits{9} = 'Tobacco packs/year';

        ScreeningData.alcohol = data(:,8);
        ScreeningData.Properties.VariableDescriptions{10} = 'Drinks alcohol? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.alcoholeh = data(:,9);
        ScreeningData.Properties.VariableUnits{11} = 'EtOH bottles/day';

        ScreeningData.hypert = data(:,10);
        ScreeningData.Properties.VariableDescriptions{12} = 'Hypertension? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.hypmed = data(:,11);
        ScreeningData.Properties.VariableDescriptions{13} = 'Medication for hypertension? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.hypmedn = data(:,12);
        ScreeningData.Properties.VariableUnits{11} = 'Units of medication';
        ScreeningData.Properties.VariableDescriptions{14} = 'How many types of medications for hypertension?';

        ScreeningData.diabet = data(:,13);
        ScreeningData.Properties.VariableDescriptions{15} = 'Diabetes? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.diahba1c = data(:,14);
        ScreeningData.Properties.VariableUnits{16} = 'mmol HbA1c/mol';
        ScreeningData.Properties.VariableDescriptions{16} = 'glycated haemoglobin (A1c) concentration';

        ScreeningData.diamed = data(:,15);
        ScreeningData.Properties.VariableDescriptions{17} = 'Medication for diabetes? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.oraal = data(:,16);
        ScreeningData.Properties.VariableDescriptions{18} = 'Oral anti-diabetic medication 0=No, 1=Yes, 9=Unknown';

        ScreeningData.oraaln = data(:,17);
        ScreeningData.Properties.VariableUnits{19} = 'Units of medication';
        ScreeningData.Properties.VariableDescriptions{19} = 'How many types of anti-diabetics?';

        ScreeningData.insul = data(:,18);
        ScreeningData.Properties.VariableDescriptions{20} = 'Insulin dependent? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.insuln = data(:,19);
        ScreeningData.Properties.VariableUnits{21} = 'Species of medication';
        ScreeningData.Properties.VariableDescriptions{21} = 'How many species (insulin) medicines?';

        ScreeningData.insuleh = data(:,20);
        ScreeningData.Properties.VariableUnits{22} = 'Units of medication';
        ScreeningData.Properties.VariableDescriptions{22} = 'How many (insulin) units total?';

        ScreeningData.dyslip = data(:,21);
        ScreeningData.Properties.VariableDescriptions{23} = 'Dyslipidemia? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.dyslipldl = data(:,22);
        ScreeningData.Properties.VariableUnits{24} = 'mmol/L';
        ScreeningData.Properties.VariableDescriptions{24} = 'Low Density Lipoprotein (LDL) concentration';

        ScreeningData.dysliphdl = data(:,23);
        ScreeningData.Properties.VariableUnits{25} = 'mmol/L';
        ScreeningData.Properties.VariableDescriptions{25} = 'High Density Lipoprotein (HDL) concentration';

        ScreeningData.dysliptrig = data(:,24);
        ScreeningData.Properties.VariableUnits{26} = 'mmol/L';
        ScreeningData.Properties.VariableDescriptions{26} = 'Triglycerides concentration';

        ScreeningData.dyslipratio = data(:,25);
        ScreeningData.Properties.VariableDescriptions{27} = 'LDL/HDL';

        ScreeningData.dyslipmed = data(:,26);
        ScreeningData.Properties.VariableDescriptions{28} = 'Medications for dyslipidemia? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.dyslipmedn = data(:,27);
        ScreeningData.Properties.VariableUnits{29} = 'Units of medication';
        ScreeningData.Properties.VariableDescriptions{29} = 'How many types of medications for dyslipidemia?';

        ScreeningData.oesofa = data(:,28);
        ScreeningData.Properties.VariableDescriptions{30} = 'Gastro-esophageal-reflux-disease (GERD) 0=No, 1=Yes, 9=Unknown';

        ScreeningData.oesofadiag = data(:,29);
        ScreeningData.Properties.VariableDescriptions{31} = 'How is it confirmed? 1=gastroduodenoscopy, 2=pH measurement, 9=unknown';

        ScreeningData.oesofamed = data(:,30);
        ScreeningData.Properties.VariableDescriptions{32} = 'Medication for GERD? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.osas = data(:,31);
        ScreeningData.Properties.VariableDescriptions{33} = 'Sleep apnea or obesity-related hypoventilation (OSAS) 0=No, 1=Yes, 9=Unknown';

        ScreeningData.cpap = data(:,32);
        ScreeningData.Properties.VariableDescriptions{34} = 'Continuous Positive Airway Pressure (CPAP) 0=No, 1=Yes, 9=Unknown';

        ScreeningData.gewrklacht = data(:,33);
        ScreeningData.Properties.VariableDescriptions{35} = 'Arthralgia? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comorb = data(:,34);
        ScreeningData.Properties.VariableDescriptions{36} = 'Is there other comorbidities? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comorbcar = data(:,35);
        ScreeningData.Properties.VariableDescriptions{37} = 'Are there any cardiac abnormalities? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comorbvas = data(:,36);
        ScreeningData.Properties.VariableDescriptions{38} = 'Are there vascular abnormalities? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comorbpul = data(:,37);
        ScreeningData.Properties.VariableDescriptions{39} = 'Are there pulmonar abnormalities? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.commda = data(:,38);
        ScreeningData.Properties.VariableDescriptions{40} = 'Are there esophagus/stomach disorders? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comhep = data(:,39);
        ScreeningData.Properties.VariableDescriptions{41} = 'Are there hepatobiliary disorders? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.psych = data(:,40);
        ScreeningData.Properties.VariableDescriptions{42} = 'Are there psychiatric disorders? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comtro1 = data(:,41);
        ScreeningData.Properties.VariableDescriptions{43} = 'Is there DVT/Pulmonary embolism? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comcar04 = data(:,42);
        ScreeningData.Properties.VariableDescriptions{44} = 'heart valve defects? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comcar01 = data(:,43);
        ScreeningData.Properties.VariableDescriptions{45} = 'Myocardial infarction? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comcar08 = data(:,44);
        ScreeningData.Properties.VariableDescriptions{46} = 'Congestive heart failure? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.corafw = data(:,45);
        ScreeningData.Properties.VariableDescriptions{47} = 'Coronary deviation? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comcar07 = data(:,46);
        ScreeningData.Properties.VariableDescriptions{48} = 'Arrhythmia? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.compul1 = data(:,47);
        ScreeningData.Properties.VariableDescriptions{49} = 'Chronic obstructive pulmonary disease (COPD)? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.fibro = data(:,48);
        ScreeningData.Properties.VariableDescriptions{50} = 'Pulmonary Fibrosis? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.compul3 = data(:,49);
        ScreeningData.Properties.VariableDescriptions{51} = 'Lung resection? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.schiz = data(:,50);
        ScreeningData.Properties.VariableDescriptions{52} = 'Schizophrenia? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.depri = data(:,51);
        ScreeningData.Properties.VariableDescriptions{53} = 'Depression? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.bipol = data(:,52);
        ScreeningData.Properties.VariableDescriptions{54} = 'Bipolar disorder? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.angst = data(:,53);
        ScreeningData.Properties.VariableDescriptions{55} = 'Anxiety Disorder? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.menstru = data(:,54);
        ScreeningData.Properties.VariableDescriptions{56} = 'Are there menstrual abnormalities? 0=No, 1=Yes, 9=Unknown';

        ScreeningData.comorbove = data(:,55);
        ScreeningData.Properties.VariableDescriptions{57} = 'Other relevant comorbidity ';

        ScreeningData.laparovg = data(:,56);
        ScreeningData.Properties.VariableDescriptions{58} = 'Laparotomy in history?';

        ScreeningData.days_from_ok = data(:,57);
        ScreeningData.Properties.VariableUnits{59} = 'days';
        ScreeningData.Properties.VariableDescriptions{59} = 'Surgery date - datbez1';

    % Clear temporary variables
        clearvars data raw raw0_0 raw0_1 raw0_2 raw0_3 raw0_4 cellVectors R ...
            dates dates0_0 dates0_1 dates0_2 dates0_3 dates0_4;

end

% Replace "999" by NaN in buikomv and hooggew (999=unknown)
    buikomv_999 = ScreeningData.buikomv==999;
    buikomv_0 = ScreeningData.buikomv==0;
    ScreeningData.buikomv(buikomv_999)= NaN;
    ScreeningData.buikomv(buikomv_0)= NaN;
    
    hooggew_999 = ScreeningData.hooggew==999;
    ScreeningData.hooggew(hooggew_999)= NaN;

% Replace 9 by NaN for other variables (9=unknown)
    ff = {'alcohol','hypert','hypmed','diabet','diamed','oraal',...
        'insul','dyslip','dyslipmed','oesofa','oesofadiag','oesofamed',...
        'osas','cpap','gewrklacht','comorb','comorbcar','comorbvas',...
        'comorbpul','commda','comhep','psych','comtro1','comcar04',...
        'comcar01','comcar08','corafw','comcar07','compul1','fibro',...
        'compul3','schiz','depri','bipol','angst','menstru',...
        'comorbove','laparovg'};
    
    for pp=1:size(ff,2)
        var_name = ff{pp};
        idx_9 = ismissing(ScreeningData(:,var_name),9);
        ScreeningData(idx_9,var_name) = array2table(NaN);
    end
    clearvars pp var_name idx_9 
    
% Replace NaN's by zero from variable 8 "roken"
    a=table2array(ScreeningData(:,8:end));
    a(isnan(a))=0;
    ScreeningData(:,8:end)=array2table(a);

% Delete variables performed in Lab Analyses
    % Values stored in the LAB file is the main source.
    % Values are copied to e.g. DATO (manually).
    ScreeningData.dyslipldl=[]; % LDL (mmol/L)
    ScreeningData.dysliphdl=[]; % HDL (mmol/l)
    ScreeningData.dysliptrig=[]; % Triglycerides (mmol/L)
    ScreeningData.diahba1c=[]; % HbA1c (mmol/L)
    % Values estimated based on the previous entries
    ScreeningData.dyslipratio=[]; % LDL/HDL- ratio

% Save raw data
    save([subfold,'/',d,'/RawData/ScreeningDataRAW'],'ScreeningData');
    
% Remove patients with surgeries before 01-01-2012
    ScreeningData(~ismember([ScreeningData.PatientCode], total_patientsID),:)=[];

clearvars a buikomv_0 buikomv_999 hooggew_999 i ff opts;
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%    Workbook: DATO CZE 2016-03.xlsx
%    Worksheet: verrichting (Operation)

if TF == 1
    % Setup the Import Options
        opts = spreadsheetImportOptions("NumVariables", 55);

    % Specify sheet and range
        opts.Sheet = "verrichting";
        opts.DataRange = "A2:BC3620";

    % Specify column names and types
        opts.VariableNames = ["PatientCode", "datok", "operateur", "operateur2", "asascore", "ingreep", "procok", "techother", "bandtype", "bandtypother", "bandtech", "bandfix", "sleebou", "sleetech", "andersnl", "afmtrans", "byintest", "byloopbili", "byloopali", "petclosed", "jejuclosed", "loopali", "loopcom", "benad", "cali", "filling", "peropcomp", "peropcomlsetall", "comper", "combloe", "commilt", "comlev", "comvat", "comdood", "compl", "compchir", "aardcomchir", "aardcomchirnl", "algcompl", "aardcomalg", "aardcomalgnl", "icopname", "redic", "reintreq", "reinttyp", "anestreint", "letsel", "datont", "heropn", "datheropn", "redheropn", "datheropnont", "status", "datovl", "doodoorz"];
        opts.SelectedVariableNames = ["PatientCode", "datok", "operateur", "operateur2", "asascore", "ingreep", "procok", "techother", "bandtype", "bandtypother", "bandtech", "bandfix", "sleebou", "sleetech", "andersnl", "afmtrans", "byintest", "byloopbili", "byloopali", "petclosed", "jejuclosed", "loopali", "loopcom", "benad", "cali", "filling", "peropcomp", "peropcomlsetall", "comper", "combloe", "commilt", "comlev", "comvat", "comdood", "compl", "compchir", "aardcomchir", "aardcomchirnl", "algcompl", "aardcomalg", "aardcomalgnl", "icopname", "redic", "reintreq", "reinttyp", "anestreint", "letsel", "datont", "heropn", "datheropn", "redheropn", "datheropnont", "status", "datovl", "doodoorz"];
        opts.VariableTypes = ["double", "datetime", "double", "string", "double", "double", "double", "string", "string", "string", "string", "string", "double", "double", "string", "double", "string", "string", "string", "string", "string", "string", "string", "double", "string", "string", "double", "string", "string", "string", "string", "string", "string", "string", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "datetime", "string", "string", "string", "string", "string", "string", "string"];
        opts = setvaropts(opts, 2, "InputFormat", "");
        opts = setvaropts(opts, 48, "InputFormat", "");
        opts = setvaropts(opts, [4, 8, 9, 10, 11, 12, 15, 17, 18, 19, 20, 21, 22, 23, 25, 26, 28, 29, 30, 31, 32, 33, 34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53, 54, 55], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, [4, 8, 9, 10, 11, 12, 15, 17, 18, 19, 20, 21, 22, 23, 25, 26, 28, 29, 30, 31, 32, 33, 34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53, 54, 55], "EmptyFieldRule", "auto");

    % Import the data
        SurgeryData = readtable("/Volumes/GoogleDrive/My Drive/MIT-Portugal/PhD Work/Bariatrics/BariatricsModel/ExtractionScripts/DATO CZE 2016-03.xlsx", opts, "UseExcel", false);

    clear opts
    
    % Allocate imported array to column variable names
        SurgeryData.Properties.VariableDescriptions{1} = 'Patient ID';
        SurgeryData.Properties.VariableDescriptions{2} = 'Date of surgery';
        SurgeryData.Properties.VariableDescriptions{3} = 'Surgeon';
        SurgeryData.Properties.VariableDescriptions{4} = 'ASA class (I, II, III, IV, V)';

else
    % Import data from spreadsheet
        [~, ~, raw0_0, dates0_0] = xlsread('DATO CZE 2016-03.xlsx','verrichting','A2:C3620','',@convertSpreadsheetExcelDates);
        [~, ~, raw0_1, dates0_1] = xlsread('DATO CZE 2016-03.xlsx','verrichting','E2:E3620','',@convertSpreadsheetExcelDates);
    raw = [raw0_0,raw0_1];
    dates = [dates0_0,dates0_1];
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    raw = raw(:,[1,3,4]);
    dates = dates(:,2);

    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
    dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN

    % Create output variable
    data = reshape([raw{:}],size(raw));

    % Create table
        SurgeryData = table;

    % Allocate imported array to column variable names
        SurgeryData.PatientCode = data(:,1);
        SurgeryData.Properties.VariableDescriptions{1} = 'Patient ID';

        SurgeryData.datok = datenum(datetime([dates{:,1}].', 'ConvertFrom', 'Excel', 'Format', 'dd/MM/yyyy'));
        SurgeryData.Properties.VariableDescriptions{2} = 'Date of surgery';

        SurgeryData.operateur = data(:,2);
        SurgeryData.Properties.VariableDescriptions{3} = 'Surgeon';

        SurgeryData.asascore = data(:,3);
        SurgeryData.Properties.VariableDescriptions{4} = 'ASA class (I, II, III, IV, V)';

    % Clear temporary variables
    clearvars data raw dates raw0_0 dates0_0 raw0_1 dates0_1 cellVectors R;
end

% Save raw data
    save([subfold,'/',d,'/RawData/SurgeryDataRAW'],'SurgeryData');

% Remove patients with surgeries before 01-01-2012
    SurgeryData(~ismember([SurgeryData.PatientCode], total_patientsID),:)=[];

%% Import follow-up data, extracting spreadsheet dates in Excel serial date format
% Script for importing data from the following spreadsheet:
%    Original Workbook: DATO CZE 2016-03.xlsx
%    Worksheet: fup
%    Extracted from: Follow-upDATO CZE 2017.xlsx
%    Worksheet: fup (contains a colum with days_from_ok)

if TF == 1
    % Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 68);

    % Specify sheet and range
    opts.Sheet = "fup";
    opts.DataRange = "A2:BP14923";

    % Specify column names and types
    opts.VariableNames = ["PatientCode", "Var2", "Var3", "followdat", "folgewicht", "hypertbeter", "Var7", "Var8", "diabbeter", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "dyslipbeter", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "osasbeter", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43", "Var44", "Var45", "Var46", "Var47", "Var48", "Var49", "Var50", "Var51", "Var52", "Var53", "Var54", "Var55", "Var56", "Var57", "Var58", "Var59", "Var60", "Var61", "Var62", "Var63", "Var64", "Var65", "Var66", "Var67", "days_from_ok"];
    opts.SelectedVariableNames = ["PatientCode", "followdat", "folgewicht", "hypertbeter", "diabbeter", "dyslipbeter", "osasbeter", "days_from_ok"];
    opts.VariableTypes = ["double", "char", "char", "datetime", "double", "double", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "string"];
    opts = setvaropts(opts, 4, "InputFormat", "");
    opts = setvaropts(opts, [2, 3, 7, 8, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [2, 3, 7, 8, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68], "EmptyFieldRule", "auto");

    % Import the data
    FollowUpData = readtable("Follow-upDATO CZE 2017.xlsx", opts, "UseExcel", false);
    
    % Allocate imported array to column variable names
        FollowUpData.Properties.VariableDescriptions{'PatientCode'} = 'Patient ID';
        FollowUpData.Properties.VariableDescriptions{'followdat'} = 'Date of Follow-up consultation';
        FollowUpData.Properties.VariableUnits{'folgewicht'} = 'kg';
        FollowUpData.Properties.VariableDescriptions{'folgewicht'} = 'Weight on Follow-up moment';
        FollowUpData.Properties.VariableDescriptions{'hypertbeter'} = 'Hypertense patient (boolean)';
        FollowUpData.Properties.VariableDescriptions{'diabbeter'} = 'Diabetic patient (boolean)';
        FollowUpData.Properties.VariableDescriptions{'dyslipbeter'} = 'Dyslipidemic patient (boolean)';
        FollowUpData.Properties.VariableDescriptions{'osasbeter'} = 'OSAS Scores: Healed/Genezen - Better/Beter - Same/Gelijk - Worse/Slechter - NA/NVT'; 
        FollowUpData.Properties.VariableUnits{'days_from_ok'} = 'days';
        FollowUpData.Properties.VariableDescriptions{'days_from_ok'} = 'Surgery date - followdat'; 

    % Clear temporary variables
        clear opts
else
    % Import data
        [~, ~, raw0_0, dates0_0] = xlsread('Follow-upDATO CZE 2017.xlsx','fup','A2:A14923','',@convertSpreadsheetExcelDates);
        [~, ~, raw0_1, dates0_1] = xlsread('Follow-upDATO CZE 2017.xlsx','fup','D2:F14923','',@convertSpreadsheetExcelDates);
        [~, ~, raw0_2, dates0_2] = xlsread('Follow-upDATO CZE 2017.xlsx','fup','I2:I14923','',@convertSpreadsheetExcelDates);
        [~, ~, raw0_3, dates0_3] = xlsread('Follow-upDATO CZE 2017.xlsx','fup','Q2:Q14923','',@convertSpreadsheetExcelDates);
        [~, ~, raw0_4, dates0_4] = xlsread('Follow-upDATO CZE 2017.xlsx','fup','X2:X14923','',@convertSpreadsheetExcelDates);
        [~, ~, raw0_5, dates0_5] = xlsread('Follow-upDATO CZE 2017.xlsx','fup','BP2:BP14923','',@convertSpreadsheetExcelDates);
    raw = [raw0_0,raw0_1,raw0_2,raw0_3,raw0_4,raw0_5];
    dates = [dates0_0,dates0_1,dates0_2,dates0_3,dates0_4,dates0_5];
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,1);
    raw = raw(:,[3,4,5,6,7,8]);
    dates = dates(:,2);

    % Replace non-numeric cells with NaN
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
        raw(R) = {NaN}; % Replace non-numeric cells
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
        dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN

    % Create output variable
        data = reshape([raw{:}],size(raw));

    % Create table
        FollowUpData = table;

    % Allocate imported array to column variable names
        FollowUpData.PatientCode = str2double(cellVectors(:,1));
        FollowUpData.Properties.VariableDescriptions{1} = 'Patient ID';

        FollowUpData.followdat = datenum(datetime([dates{:,1}].', 'ConvertFrom', 'Excel'));
        FollowUpData.Properties.VariableDescriptions{2} = 'Date of Follow-up consultation';

        FollowUpData.folgewicht = data(:,1);
        FollowUpData.Properties.VariableUnits{3} = 'kg';
        FollowUpData.Properties.VariableDescriptions{3} = 'Weight on Follow-up moment';

        FollowUpData.hypertbeter = data(:,2);
        FollowUpData.Properties.VariableDescriptions{4} = 'Hypertense patient (boolean)';

        FollowUpData.diabbeter = data(:,3);
        FollowUpData.Properties.VariableDescriptions{5} = 'Diabetic patient (boolean)';

        FollowUpData.dyslipbeter = data(:,4);
        FollowUpData.Properties.VariableDescriptions{6} = 'Dyslipidemic patient (boolean)';

        FollowUpData.osasbeter = data(:,5);
        FollowUpData.Properties.VariableDescriptions{7} = 'OSAS Scores: Healed/Genezen - Better/Beter - Same/Gelijk - Worse/Slechter - NA/NVT'; 

        FollowUpData.days_from_ok = data(:,6);
        FollowUpData.Properties.VariableUnits{8} = 'days';
        FollowUpData.Properties.VariableDescriptions{8} = 'Surgery date - followdat'; 

    % Clear temporary variables
        clearvars data raw dates raw0_0 dates0_0 raw0_1 dates0_1 raw0_2 dates0_2 ...
        raw0_3 dates0_3 raw0_4 dates0_4 raw0_5 dates0_5 cellVectors R;
end

% Save raw data
    save([subfold,'/',d,'/RawData/FollowUpDataRAW'],'FollowUpData');

% Remove patients with surgeries before 01-01-2012
FollowUpData(~ismember([FollowUpData.PatientCode], total_patientsID),:)=[];

%% Import data from spreadsheet
    % Script for importing data from the following spreadsheet:
        % Workbook: LAB 2016-03.xlsx
        % Worksheet: LabDATA
        % Deleted data before surgery (>=0)

if TF == 1
    % Setup the Import Options
        opts = spreadsheetImportOptions("NumVariables", 17);

    % Specify sheet and range
        opts.Sheet = "LabDATA";
        opts.DataRange = "A2:Q290965";

    % Specify column names and types
        opts.VariableNames = ["PatientCode", "Var2", "bepcode", "description", "uitslag", "uitslagNUM", "eenheid", "Var8", "datumtijd", "refwaarde", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "days_from_ok"];
        opts.SelectedVariableNames = ["PatientCode", "bepcode", "description", "uitslag", "uitslagNUM", "eenheid", "datumtijd", "refwaarde", "days_from_ok"];
        opts.VariableTypes = ["double", "char", "string", "string", "double", "double", "string", "char", "datetime", "string", "char", "char", "char", "char", "char", "char", "double"];
        opts = setvaropts(opts, 9, "InputFormat", "");
        opts = setvaropts(opts, [2, 3, 4, 7, 8, 10, 11, 12, 13, 14, 15, 16], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, [2, 3, 4, 7, 8, 10, 11, 12, 13, 14, 15, 16], "EmptyFieldRule", "auto");

    % Import the data
        LabResults = readtable("LAB_CZE_2017_CLEAN.xlsx", opts, "UseExcel", false);
    
    % Allocate imported array to column variable names
        LabResults.Properties.VariableDescriptions{'PatientCode'} = 'Patient ID';

        LabResults.Properties.VariableDescriptions{'bepcode'} = 'Analysis ID';

        LabResults.Properties.VariableDescriptions{'description'} = 'Analysis name';

        LabResults.Properties.VariableDescriptions{'uitslag'} = 'Result';

        LabResults.Properties.VariableDescriptions{'uitslagNUM'} = 'Result';

        LabResults.Properties.VariableDescriptions{'datumtijd'} = 'Date of analysis';

        LabResults.Properties.VariableDescriptions{'refwaarde'} = 'Reference value(s)';

        LabResults.Properties.VariableUnits{'days_from_ok'} = 'days';
        LabResults.Properties.VariableDescriptions{'days_from_ok'} = 'measurement date - surgery date (numerical)';
    
    % Clear temporary variables
        clear opts
    
else
    % Import the data, extracting spreadsheet dates in Excel serial date format
        [~, ~, raw0_0, dates0_0] = xlsread('LAB_CZE_2017_CLEAN.xlsx','LabDATA','A2:A290965','',@convertSpreadsheetExcelDates);
        [~, ~, raw0_1, dates0_1] = xlsread('LAB_CZE_2017_CLEAN.xlsx','LabDATA','C2:G290965','',@convertSpreadsheetExcelDates);
        [~, ~, raw0_2, dates0_2] = xlsread('LAB_CZE_2017_CLEAN.xlsx','LabDATA','I2:J290965','',@convertSpreadsheetExcelDates);
        [~, ~, raw0_3, dates0_3] = xlsread('LAB_CZE_2017_CLEAN.xlsx','LabDATA','Q2:Q290965','',@convertSpreadsheetExcelDates);
    raw = [raw0_0,raw0_1,raw0_2,raw0_3];
    dates = [dates0_0,dates0_1,dates0_2,dates0_3];
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,[1,2,3,4,6,8]);
    raw = raw(:,[5,9]);
    dates = dates(:,7);

    % Replace non-numeric cells with NaN
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
        raw(R) = {NaN}; % Replace non-numeric cells
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
        dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN

    % Create output variable
        data = reshape([raw{:}],size(raw));

    % Create table
        LabResults = table;

    % Allocate imported array to column variable names
        LabResults.PatientCode = str2double(cellVectors(:,1));
        LabResults.Properties.VariableDescriptions{1} = 'Patient ID';

        LabResults.bepcode = cellVectors(:,2);
        LabResults.Properties.VariableDescriptions{2} = 'Analysis ID';

        LabResults.description = cellVectors(:,3);
        LabResults.Properties.VariableDescriptions{3} = 'Analysis name';

        LabResults.uitslag = cellVectors(:,4);
        LabResults.Properties.VariableDescriptions{4} = 'Result';

        LabResults.uitslagNUM = data(:,1);
        LabResults.Properties.VariableDescriptions{5} = 'Result';

        LabResults.units = cellVectors(:,5);

        LabResults.datumtijd = datenum(datetime([dates{:,1}].', 'ConvertFrom', 'Excel', 'Format', 'dd/MM/yyyy'));
        LabResults.Properties.VariableDescriptions{7} = 'Date of analysis';

        LabResults.refwaarde = cellVectors(:,6);
        LabResults.Properties.VariableDescriptions{8} = 'Reference value(s)';

        LabResults.days_from_ok = data(:,2);
        LabResults.Properties.VariableUnits{9} = 'days';
        LabResults.Properties.VariableDescriptions{9} = 'measurement date - surgery date (numerical)';

    % Clear temporary variables
        clearvars data raw dates raw0_0 dates0_0 raw0_1 dates0_1 raw0_2 dates0_2 raw0_3 dates0_3 cellVectors R;
end

% Removed data
% Data before surgery days_from_ok<0
% Data with bepcode count less than 100
% Remove data containing info about blood type
%     LabResults(~cellfun(@isempty,strfind(table2array(LabResults(:,'bepcode')), 'BBL')),:)=[];
% 
% Remove data collected with old lab equipment
%     LabResults(~cellfun(@isempty,strfind(table2array(LabResults(:,'description')), '(oud)')),:)=[];
% 
% Remove codes that contain no lab results
%     LabResults(~cellfun(@isempty,strfind(table2array(LabResults(:,'bepcode')), 'REM')),:)=[];
% 
% Remove other irrelevant lab results
%     LabResults(strcmp(LabResults.bepcode, 'SCO000'),:)=[];  % Cortisol in speeksel
%     LabResults(strcmp(LabResults.bepcode, 'BIR000'),:)=[];  % Irr. antistoffen screening
%     LabResults(strcmp(LabResults.bepcode, 'BME010'),:)=[];  % Medicijnspiegel
%     LabResults(strcmp(LabResults.bepcode, 'BSP000'),:)=[];  % Spijtserum
%     LabResults(strcmp(LabResults.bepcode, 'REM3'),:)=[];    % Telefonisch melden aan
%     LabResults(strcmp(LabResults.bepcode, 'XTR003'),:)=[];  % Trf rapportage
%     LabResults(strcmp(LabResults.bepcode, 'BVA004'),:)=[];  % Vancomycine dal
%     LabResults(strcmp(LabResults.bepcode, 'BVA003'),:)=[];  % Vancomycine top

LabResults(strcmp(LabResults.bepcode, 'TGL001'),:)=[];
LabResults(strcmp(LabResults.bepcode, 'ULE004'),:)=[];
% LabResults(strcmp(LabResults.bepcode, 'BVI032'),:)=[]; % VitA 
LabResults(strcmp(LabResults.bepcode, 'XCO000'),:)=[];
LabResults(strcmp(LabResults.bepcode, 'USC002'),:)=[];
% LabResults(strcmp(LabResults.bepcode, 'BVI004'),:)=[]; % VitB6_BV
LabResults(strcmp(LabResults.bepcode, 'BVI015'),:)=[]; % VitK1_BV
LabResults(strcmp(LabResults.bepcode, 'BVI021'),:)=[]; % 1_25diOH_VitD3
LabResults(strcmp(LabResults.bepcode, 'BVI024'),:)=[]; % VitB2_BV
LabResults(strcmp(LabResults.bepcode, 'BVI025'),:)=[]; % VitC_BV
LabResults(strcmp(LabResults.bepcode, 'BVI026'),:)=[]; % VitE_BV
% LabResults(strcmp(LabResults.bepcode, 'BVI030'),:)=[]; % 25OHD3VitamineD3_BV
% LabResults(strcmp(LabResults.bepcode, 'BVI031'),:)=[]; % VitB1_BV
LabResults(strcmp(LabResults.bepcode, 'BVI033'),:)=[]; % VitE_BV_C
LabResults(strcmp(LabResults.bepcode, 'BVI035'),:)=[]; % VitB2_BV_C

% Remove time biased variables
% LabResults(strcmp(LabResults.bepcode, 'BZI003'),:)=[]; % Zinc
% LabResults(strcmp(LabResults.bepcode, 'BVI028'),:)=[]; % VitB12
LabResults(strcmp(LabResults.bepcode, 'BFT004'),:)=[]; % FT4
LabResults(strcmp(LabResults.bepcode, 'BME009'),:)=[]; % Methyl malonic acid

% Convert Prothrombine time (PT) into INR
PTnormal = 13; % Laboratory's geometric mean value for normal patients (seconds)
ISI = 0.935; % International Sensitivity Index

for i=1:height(LabResults)
    BEPcode = LabResults.bepcode(i);
    PTmeas = LabResults.uitslagNUM(i);
    
    if strcmp(BEPcode, 'BPR012') == 1
        LabResults.uitslagNUM(i) = (PTmeas/PTnormal)^ISI;
    end
    
    clearvars BEPcode PTmeas;
end

clearvars ISI i PTnormal;

% Save raw data
    save([subfold,'/',d,'/RawData/LabDataRAW'],'LabResults');

% Remove patients with surgeries before 01-01-2012
    LabResults(~ismember([LabResults.PatientCode], total_patientsID),:)=[];

%% Import  Physical Measurements extracting spreadsheet dates in Excel serial date format
% Script for importing data from the following spreadsheet:
    % Original Workbook: PhysicalMeas CZE 2016-03.xlsx
    % FINAL FILE: PhysucalMeas CZE 2017.xlsx. This file contains the
    % difference between datumtijd - datok (days_from_ok column)
    % Deleted data after surgery (>= 0)

if TF == 1
    % Setup the Import Options
        opts = spreadsheetImportOptions("NumVariables", 14);

    % Specify sheet and range
        opts.Sheet = "DATA";
        opts.DataRange = "A2:N109691";

    % Specify column names and types
        opts.VariableNames = ["PatientCode", "datumtijd", "vraagID", "STELLING", "korteoms", "antwoord", "omschrijving01", "Var8", "Var9", "Var10", "Var11", "value01", "value02", "days_from_ok"];
        opts.SelectedVariableNames = ["PatientCode", "datumtijd", "vraagID", "STELLING", "korteoms", "antwoord", "omschrijving01", "value01", "value02", "days_from_ok"];
        opts.VariableTypes = ["double", "datetime", "double", "string", "string", "string", "string", "char", "char", "char", "char", "double", "double", "double"];
        opts = setvaropts(opts, 2, "InputFormat", "");
        opts = setvaropts(opts, [4, 5, 6, 7, 8, 9, 10, 11], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, [4, 5, 6, 7, 8, 9, 10, 11], "EmptyFieldRule", "auto");

    % Import the data
        PhMsData = readtable("PhysicalMeas_CZE_2017.xlsx", opts, "UseExcel", false);

    % Clear temporary variables
        clear opts
    
    % Allocate imported array to column variable names
        PhMsData.Properties.VariableDescriptions{'PatientCode'} = 'Patient ID';

        PhMsData.Properties.VariableDescriptions{'datumtijd'} = 'Date of consultation/measurement (numeric)';

        PhMsData.Properties.VariableDescriptions{'vraagID'} = '(questionID) unique number for each question';

        PhMsData.Properties.VariableDescriptions{'STELLING'} = 'Question (text)';

        PhMsData.Properties.VariableDescriptions{'antwoord'} = 'Answer (text or numerical)';

        PhMsData.Properties.VariableDescriptions{'omschrijving01'} = 'Description (text or numerical)';

        PhMsData.Properties.VariableDescriptions{'value01'} = 'Response/Result';

        PhMsData.Properties.VariableDescriptions{'value02'} = 'Response/Result';

        PhMsData.Properties.VariableUnits{'days_from_ok'} = 'days';
        PhMsData.Properties.VariableDescriptions{'days_from_ok'} = 'measurement date - surgery date (numerical)';

else
    % Import data
    [~, ~, raw0_0, dates0_0] = xlsread('PhysicalMeas_CZE_2017.xlsx','DATA','A2:G109691','',@convertSpreadsheetExcelDates);
    [~, ~, raw0_1, dates0_1] = xlsread('PhysicalMeas_CZE_2017.xlsx','DATA','L2:N109691','',@convertSpreadsheetExcelDates);
        raw = [raw0_0,raw0_1];
        dates = [dates0_0,dates0_1];
        raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
        cellVectors = raw(:,[1,4,5]);
        raw = raw(:,[3,6,7,8,9,10]);
        dates = dates(:,2);

    % Replace non-numeric cells with NaN
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
        raw(R) = {NaN}; % Replace non-numeric cells
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
        dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN

    % Create output variable
        data = reshape([raw{:}],size(raw));

    % Create table
        PhMsData = table;

    % Allocate imported array to column variable names
        PhMsData.PatientCode = str2double(cellVectors(:,1));
        PhMsData.Properties.VariableDescriptions{1} = 'Patient ID';

        PhMsData.datumtijd = datenum(datetime([dates{:,1}].', 'ConvertFrom', 'Excel', 'Format', 'dd/MM/yy HH:mm:ss'));
        PhMsData.Properties.VariableDescriptions{2} = 'Date of consultation/measurement (numeric)';

        PhMsData.vraagID = data(:,1);
        PhMsData.Properties.VariableDescriptions{3} = '(questionID) unique number for each question';

        PhMsData.STELLING = cellVectors(:,2);
        PhMsData.Properties.VariableDescriptions{4} = 'Question (text)';

        PhMsData.antwoord = data(:,2);
        PhMsData.Properties.VariableDescriptions{5} = 'Answer (text or numerical)';

        PhMsData.omschrijving01 = data(:,3);
        PhMsData.Properties.VariableDescriptions{6} = 'Description (text or numerical)';

        PhMsData.value01 = data(:,4);
        PhMsData.Properties.VariableDescriptions{7} = 'Response/Result';

        PhMsData.value02 = data(:,5);
        PhMsData.Properties.VariableDescriptions{8} = 'Response/Result';

        PhMsData.days_from_ok = data(:,6);
        PhMsData.Properties.VariableUnits{9} = 'days';
        PhMsData.Properties.VariableDescriptions{9} = 'measurement date - surgery date (numerical)';

    % Clear temporary variables
        clearvars data raw dates raw0_0 dates0_0 raw0_1 dates0_1 cellVectors R;
end

% Remove question IDs
    % Contains different text responses/results
    PhMsData(ismissing(PhMsData.vraagID, 120484),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 118497),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 129580),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 146949),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 146950),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 147999),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 118502),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 120484),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 120695),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 126722),:)=[];
    PhMsData(ismissing(PhMsData.vraagID, 129581),:)=[];
    
    % Contains different boolean values at the same time for just one question ID
    PhMsData(ismissing(PhMsData.vraagID, 118511),:)=[];
    % PhMsData(ismissing(PhMsData.vraagID, 120196),:)=[];
    
    % Contains information saved in other variables
        % BMI
        PhMsData(ismissing(PhMsData.vraagID, 118506),:)=[];
        
        % Weight
        PhMsData(ismissing(PhMsData.vraagID, 112208),:)=[]; % Few incidences
        PhMsData(ismissing(PhMsData.vraagID, 132960),:)=[]; % Few incidences
    
    % Remove variables suggested by health-care profesionals
        % MIP Maximal Inspiratory?ContinuousPressure (predicted)
        PhMsData(ismissing(PhMsData.vraagID, 130472),:)=[];
    
% Add days between smoking-quit date and measuring date
    % Import the data vraagID = 120196
    [~, ~, raw] = xlsread('PhysicalMeas_NonSmokingDays_2017.xlsx','NonSmoke','A2:C983');

    % Create output variable
    PhMsNonSmokingDays = reshape([raw{:}],size(raw));

    % Obtain indexes of data obtained after surgery
    % idx_after_surg = find(PhMsNonSmokingDays(:,3) >= 0);
    %     idx_after_surg = bsxfun(@ge,PhMsNonSmokingDays(:,3),0);
    %     PhMsNonSmokingDays(idx_after_surg,:)=[];
    
    smoke_patients = unique(PhMsNonSmokingDays(:,1));
    
    for m = 1:length(smoke_patients)
        smoking_patient = smoke_patients(m);
        
        idx_smoking_patient = [PhMsData.PatientCode] == smoking_patient;
        idy_smoking_patient = PhMsNonSmokingDays(:,1) == smoking_patient;
        
        idx_vraagID = [PhMsData.vraagID] == 126718;
        
        idx = idx_smoking_patient & idx_vraagID;
        
            PhMsData.value01(idx) = PhMsNonSmokingDays(idy_smoking_patient,3);    
    end

% Clear temporary variables
    clearvars raw idx_after_surg idx idx_smoking_patient ...
        idy_smoking_patient idx_vraagID smoke_patients smoking_patient...
        PhMsNonSmokingDays m;

% Save raw data
    save([subfold,'/',d,'/RawData/PhMsDataRAW'],'PhMsData');
    
% Remove patients with surgeries before 01-01-2012
PhMsData(~ismember([PhMsData.PatientCode], total_patientsID),:)=[];
    
% Re-scaling BARO QoL (Quality of Life) score
    % Between november 2012 and january 2013 the CZE decided to enlarge the
    % BARO QoL score from 1-10 to 1-20. However, after that they decided to
    % go back to the latter scale (1-10). These values (1-20) will be 
    % lineraly transformed to a scale ranging from 1 to 10.
    % This transformation applies for codes: 130460 Esteem Score; 
    % 130461 Physical score; 130462 Sociaal score; 130464 Work score; 
    % 130465 Sexual score
% Done on the Excel file
%     DateString1 = '01-Nov-2012';
%     DateString2 = '31-Jan-2013';
%     formatIn = 'dd-mmm-yyyy';
%     Date1 = datenum(DateString1,formatIn);
%     Date2 = datenum(DateString2,formatIn);
%     BARO_codes = [130460;130461;130462;130464;130465];
%     idx_dates = (PhMsData.datumtijd >= Date1)&(PhMsData.datumtijd <= Date2);

%% Import PPOS (Polikliniek Pre-Operatieve Screening) data from spreadsheet
% Script for importing data from the following spreadsheet:
    % Original Workbook: PreSurgicalScreening CZE 2016-03.xlsx
    % Original Spreadsheets: "PPOS DATA deel1" and "PPOS DATA deel2"
    % Difference between deel 1 and deel 2 (in English: part 1 and part 2)
        % is a result of the data being extracted in 2 batches.
        % They where merged in postgreSQL. Overlaps/duplicates where
        % deleted. FINAL FILE: PreSurgicalScreeningFUSED CZE 2016-03.xlsx

if TF == 1
    % Setup the Import Options
        opts = delimitedTextImportOptions("NumVariables", 8);

    % Specify range and delimiter
        opts.DataLines = [2, Inf];
        opts.Delimiter = ",";

    % Specify column names and types
        opts.VariableNames = ["PatientCode", "datum", "vraagid", "STELLING", "antwoord", "Value01", "days_from_ok", "questionID"];
        opts.VariableTypes = ["double", "datetime", "double", "string", "double", "double", "double", "string"];
        opts = setvaropts(opts, 2, "InputFormat", "yyyy-MM-dd");
        opts = setvaropts(opts, [4, 8], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, 5, "TrimNonNumeric", true);
        opts = setvaropts(opts, 5, "ThousandsSeparator", ",");
        opts = setvaropts(opts, [4, 8], "EmptyFieldRule", "auto");
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
    
    % Import the data
        PPOSdata = readtable("PPOSdataRAW.csv", opts);

    % Clear temporary variables
        clear opts
    
    % Allocate imported array to column variable names
        PPOSdata.Properties.VariableDescriptions{'PatientCode'} = 'Patient ID';

        PPOSdata.Properties.VariableDescriptions{'datum'} = 'Date of consultation for screening (anesthesiologists) (numeric)';

        PPOSdata.Properties.VariableDescriptions{'questionID'} = '(questionID) is a unique number for each question, intencionally the vID prefix was added to create column titles';

        PPOSdata.Properties.VariableDescriptions{'STELLING'} = 'Question (text)';

        PPOSdata.Properties.VariableDescriptions{'antwoord'} = 'Answer (text or numerical)';

        PPOSdata.Properties.VariableDescriptions{'Value01'} = 'OBJECTID VALUE01';

        PPOSdata.Properties.VariableUnits{'days_from_ok'} = 'days';
        PPOSdata.Properties.VariableDescriptions{8} = 'Date of Surgery - datum';
        
    % Contains different text responses/results
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120837'),:)=[];

        % Delete variables repeated on Screening and Physical Measures datasets
        % Do you smoke?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID108975'),:)=[];

        % Did you have a heart attack?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120845'),:)=[];

        % Do you have a heart valve defect/disease?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120849'),:)=[];

        % Do you have asthma / COPD (chronic obstructive pulmonary disease)?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120858'),:)=[];

        % Do you have diabetes?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120863'),:)=[];

else
    % Import the data, extracting spreadsheet dates in Excel serial date format
            [~, ~, raw0_0, dates0_0] = xlsread('PreSurgicalScreeningFINAL_CZE_2017.xlsx','PreSurgicalScreeningFINAL_CZE_2','A2:B98575','',@convertSpreadsheetExcelDates);
            [~, ~, raw0_1, dates0_1] = xlsread('PreSurgicalScreeningFINAL_CZE_2017.xlsx','PreSurgicalScreeningFINAL_CZE_2','D2:D98575','',@convertSpreadsheetExcelDates);
            [~, ~, raw0_2, dates0_2] = xlsread('PreSurgicalScreeningFINAL_CZE_2017.xlsx','PreSurgicalScreeningFINAL_CZE_2','F2:H98575','',@convertSpreadsheetExcelDates);
            [~, ~, raw0_3, dates0_3] = xlsread('PreSurgicalScreeningFINAL_CZE_2017.xlsx','PreSurgicalScreeningFINAL_CZE_2','L2:M98575','',@convertSpreadsheetExcelDates);
        raw = [raw0_0,raw0_1,raw0_2,raw0_3];
        dates = [dates0_0,dates0_1,dates0_2,dates0_3];
        raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
        cellVectors = raw(:,[3,6,8]);
        raw = raw(:,[1,4,5,7]);
        dates = dates(:,2);

    % Replace non-numeric cells with NaN
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
        raw(R) = {NaN}; % Replace non-numeric cells
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
        dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN

    % Create output variable
        data = reshape([raw{:}],size(raw));

    % Create table
        PPOSdata = table;

    % Allocate imported array to column variable names
        PPOSdata.PatientCode = data(:,1);
        PPOSdata.Properties.VariableDescriptions{1} = 'Patient ID';

        PPOSdata.datum = datetime([dates{:,1}].', 'ConvertFrom', 'Excel', 'Format', 'dd/MM/yyyy');
        PPOSdata.Properties.VariableDescriptions{2} = 'Date of consultation for screening (anesthesiologists) (numeric)';

        PPOSdata.questionID = cellVectors(:,3);
        PPOSdata.Properties.VariableDescriptions{3} = '(questionID) is a unique number for each question, intencionally the vID prefix was added to create column titles';

        PPOSdata.STELLING = cellVectors(:,1);
        PPOSdata.Properties.VariableDescriptions{4} = 'Question (text)';

        PPOSdata.antwoord = data(:,2);
        PPOSdata.Properties.VariableDescriptions{5} = 'Answer (text or numerical)';

        PPOSdata.omschrijving01 = cellVectors(:,2);
        PPOSdata.Properties.VariableDescriptions{6} = 'Description (text or numerical)';

        PPOSdata.Value01 = data(:,3);
        PPOSdata.Properties.VariableDescriptions{7} = 'OBJECTID VALUE01';

        PPOSdata.days_from_ok = data(:,4);
        PPOSdata.Properties.VariableUnits{8} = 'days';
        PPOSdata.Properties.VariableDescriptions{8} = 'Date of Surgery - datum';

    % Clear temporary variables
        clearvars data raw dates raw0_0 dates0_0 raw0_1 dates0_1 raw0_2 ...
            dates0_2 raw0_3 dates0_3 cellVectors R;
    
    % Remove question IDs
        % Contains different text responses/results
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120837'),:)=[];

        % Delete variables repeated on Screening and Physical Measures datasets
        % Do you smoke?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID108975'),:)=[];

        % Did you have a heart attack?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120845'),:)=[];

        % Do you have a heart valve defect/disease?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120849'),:)=[];

        % Do you have asthma / COPD (chronic obstructive pulmonary disease)?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120858'),:)=[];

        % Do you have diabetes?
        PPOSdata(ismissing(PPOSdata.questionID, 'vID120863'),:)=[];

clearvars idx_antwoord idx_value01 idx_nan
end

% Remove unneeded data
    % Delete rows without numeric results
    idx_antwoord = ismissing(PPOSdata.antwoord,NaN);
    idx_value01 = ismissing(PPOSdata.Value01,NaN);
    %idx_oms01 = ismissing(PPOSdata.omschrijving01,'NaN');

    idx_nan = idx_antwoord & idx_value01;
    PPOSdata(idx_nan,:)=[];


% Convert numerical responses to categorical
    % Responses that imply an afirmative answer (YES = 1/true)
    idx_100004193=[PPOSdata.antwoord]==100004193;
    PPOSdata.antwoord(idx_100004193) = true;
    
    idx_100010480=[PPOSdata.antwoord]==100010480;
    PPOSdata.antwoord(idx_100010480) = true;
    
    idx_100031744=[PPOSdata.antwoord]==100031744;
    PPOSdata.antwoord(idx_100031744) = true;
    
    idx_100049587=[PPOSdata.antwoord]==100049587;
    PPOSdata.antwoord(idx_100049587) = true;
    
    clearvars idx_100004193 idx_100010480 idx_100031744 idx_100049587
    
    % Responses that imply a negative answer (NO = 0/false)
    idx_100004194=[PPOSdata.antwoord]==100004194;
    PPOSdata.antwoord(idx_100004194) = false;
    
    idx_100010479=[PPOSdata.antwoord]==100010479;
    PPOSdata.antwoord(idx_100010479) = false;
    
    idx_100031743=[PPOSdata.antwoord]==100031743;
    PPOSdata.antwoord(idx_100031743) = false;
    
    idx_100049586=[PPOSdata.antwoord]==100049586;
    PPOSdata.antwoord(idx_100049586) = false;
    
    % Responses regarding allergies qID = 121191
    % NONE = 0
    idx_100071678=[PPOSdata.antwoord]==100071678;
    PPOSdata.antwoord(idx_100071678) = 0;
    
    % Iodine = 1
    idx_100049853=[PPOSdata.antwoord]==100049853;
    PPOSdata.antwoord(idx_100049853) = 1;
    
    % Pleisters = 2
    idx_100049854=[PPOSdata.antwoord]==100049854;
    PPOSdata.antwoord(idx_100049854) = 2;
    
    % Latex	= 3
    idx_100049855=[PPOSdata.antwoord]==100049855;
    PPOSdata.antwoord(idx_100049855) = 3;
    
    % Medicines = 4
    idx_100049856=[PPOSdata.antwoord]==100049856;
    PPOSdata.antwoord(idx_100049856) = 4;
    
    % Other	= 5
    idx_100049858=[PPOSdata.antwoord]==100049858;
    PPOSdata.antwoord(idx_100049858) = 5;
    
    clearvars idx_100049853 idx_100049854 idx_100049855 idx_100049856 ...
        idx_100049858 idx_100071678 idx_100004194 idx_100010479 ...
        idx_100031743 idx_100049586
% Save raw data
    save([subfold,'/',d,'/RawData/PPOSdataRAW2'],'PPOSdata');
    
% Remove patients with surgeries before 01-01-2012
    PPOSdata(~ismember([PPOSdata.PatientCode], total_patientsID),:)=[];
%% Histogram to determine cut off points
    %c=-365:1:50; % Number of bins
    figure('Name','Pre-Operative Screening');
    %histogram(table2array(PPOSdata(:,'days_from_ok')),c)
    ksdensity(table2array(PPOSdata(:,'days_from_ok'))...
        ,'Censoring',(table2array(PPOSdata(:,'days_from_ok'))<-365));
    xL=get(gca,'Xlim');
    yL=get(gca,'Ylim');
    line([-30 -30], yL, 'Color','#A2142F','LineStyle','--','LineWidth',1.5,'DisplayName','-1 month');
    line([-60 -60], yL, 'Color','#77AC30','LineStyle','--','LineWidth',1.5,'DisplayName','-2 months');
    line([-90 -90], yL, 'Color','#7E2F8E','LineStyle','--','LineWidth',1.5,'DisplayName','-3 months');
    xlabel('Days before surgery','FontSize',16,'FontWeight','bold')
    ylabel('Density','FontSize',16,'FontWeight','bold')
    xlim([xL(:,1) 0]);
    title('Pre-Operative Screening','FontSize',18)
    legend({'Data','-1 month','-2 months','-3 months'},'Location','best','FontSize',14)
    fig = gcf; ax = gca;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 28; %Width
    fig.Position(4) = 20; %Height
    ax.FontSize = 14;
    print([subfold,'/',d,'/Plots/PreOpScreenFreq'],'-deps')
    print([subfold,'/',d,'/Plots/PreOpScreenFreqc'],'-depsc')
    print([subfold,'/',d,'/Plots/PreOpScreenFreq'], '-dpdf', '-fillpage')
    clear c xL yL

    figure('Name','Lab Analysis');
    %c=-365:1:0;
    %histogram(table2array(LabResults(:,'days_from_ok')),c)
    ksdensity(table2array(LabResults(:,'days_from_ok'))...
        ,'Censoring',((table2array(LabResults(:,'days_from_ok')))<-365)&...
        (table2array(LabResults(:,'days_from_ok'))<=0));
    yL=get(gca,'Ylim');
    line([-30 -30], yL, 'Color','#A2142F','LineStyle','--','LineWidth',1.5,'DisplayName','-1 month');
    line([-60 -60], yL, 'Color','#77AC30','LineStyle','--','LineWidth',1.5,'DisplayName','-2 months');
    line([-90 -90], yL, 'Color','#7E2F8E','LineStyle','--','LineWidth',1.5,'DisplayName','-3 months');
    xlabel('Days before surgery','FontSize',16,'FontWeight','bold')
    ylabel('Density','FontSize',16,'FontWeight','bold')
    title('Lab tests','FontSize',18)
    legend({'Data','-1 month','-2 months','-3 months'},'Location','best','FontSize',14)
    fig = gcf; ax = gca;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 28; %Width
    fig.Position(4) = 20; %Height
    ax.FontSize = 14;
    print([subfold,'/',d,'/Plots/LabFreq'],'-deps')
    print([subfold,'/',d,'/Plots/LabFreqc'],'-depsc')
    print([subfold,'/',d,'/Plots/LabFreq'], '-dpdf', '-fillpage')
    clear c yL
    
    figure('Name','Physical Measures');
    %c=-365:1:0;
    %histogram(table2array(PhMsData(:,'days_from_ok')),c)
    ksdensity(table2array(PhMsData(:,'days_from_ok'))...
        ,'Censoring',(table2array(PhMsData(:,'days_from_ok'))<-365));
    yL=get(gca,'Ylim');
    line([-30 -30], yL, 'Color','#A2142F','LineStyle','--','LineWidth',1.5,'DisplayName','-1 month');
    line([-60 -60], yL, 'Color','#77AC30','LineStyle','--','LineWidth',1.5,'DisplayName','-2 months');
    line([-90 -90], yL, 'Color','#7E2F8E','LineStyle','--','LineWidth',1.5,'DisplayName','-3 months');
    xlabel('Days before surgery','FontSize',16,'FontWeight','bold')
    ylabel('Density','FontSize',16,'FontWeight','bold')
    title('Physical Measures','FontSize',18)
    legend({'Data','-1 month','-2 months','-3 months'},'Location','best','FontSize',14)
    fig = gcf; ax = gca;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 28; %Width
    fig.Position(4) = 20; %Height
    ax.FontSize = 14;
    print([subfold,'/',d,'/Plots/PhMFreq'],'-deps')
    print([subfold,'/',d,'/Plots/PhMFreqc'],'-depsc')
    print([subfold,'/',d,'/Plots/PhMFreq'], '-dpdf', '-fillpage')
    clear c yL
    
    figure('Name','Follow-up ');
    c=0:10:((2*365)+(30*6)); % up to 2 years + 6 months
    if TF == 1
        A = table2cell(FollowUpData(:,'days_from_ok'));
        A = cellfun(@str2double,A);
        histogram(A,c,'DisplayName','Data','FaceColor','#EDB120',...
            'EdgeColor','none')
    else
        histogram(table2array(FollowUpData(:,'days_from_ok')),c,...
            'DisplayName','Data')
    end
    xL=get(gca,'Ylim');
%     line([180 180], xL, 'Color',[0 0.4470 0.7410],'LineStyle','--','DisplayName','1 year - 6 months','LineWidth',1.5); % 1 year - 6 months
    line([275 275], xL, 'Color',[0.8500 0.3250 0.0980],'LineStyle','--','DisplayName','1 year - 3 months','LineWidth',1.5); % 1 year - 3 months
    line([365 365], xL, 'Color','#A2142F','LineStyle','--','DisplayName','1 year','LineWidth',1.5); % 1 year
    line([455 455], xL, 'Color',[0.8500 0.3250 0.0980],'LineStyle','--','DisplayName','1 year + 3 months','LineWidth',1.5); % 1 year + 3 months
%     line([545 545], xL, 'Color',[0.4940 0.1840 0.5560],'LineStyle','--','DisplayName','1 year + 6 months','LineWidth',1.5); % 1 year + 6 months
    line([640 640], xL, 'Color',[0.4660 0.6740 0.1880],'LineStyle','--','DisplayName','2 years - 3 months','LineWidth',1.5); % 2 years - 3 months
    line([730 730], xL, 'Color',[0.3010 0.7450 0.9330],'LineStyle','--','DisplayName','2 years','LineWidth',1.5); % 2 years
    line([820 820], xL, 'Color','#77AC30','LineStyle','--','DisplayName','2 years + 3 months','LineWidth',1.5); % 2 years + 3 months
%     line([910 910], xL, 'Color',[0.6350 0.0780 0.1840],'LineStyle','--','DisplayName','2 years + 6 months','LineWidth',1.5); % 2 years + 6 months
    xlabel('Days after surgery','FontSize',16,'FontWeight','bold')
    ylabel('Frequency','FontSize',16,'FontWeight','bold')
    legend('Location','bestoutside','FontSize',14)
    if TF ==1 ;xlim([0 max(c)+10]);end
    fig = gcf; ax = gca;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 28; %Width
    fig.Position(4) = 20; %Height
    ax.FontSize = 14;
    print([subfold,'/',d,'/Plots/FollowUpFreq'],'-deps')
    print([subfold,'/',d,'/Plots/FollowUpFreqc'],'-depsc')
    print([subfold,'/',d,'/Plots/FollowUpFreq'], '-dpdf', '-fillpage')
    clear c xL yL TF
%% Created by Aldo Arévalo