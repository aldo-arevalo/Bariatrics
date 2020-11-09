%% Standardizing lab values using their reference intervals
% Code by Caro Fucs
% Adapted by Aldo Arevalo
%
% The code uses the structure Patient with field LabResults as input.
% The field LabResults is a table containing at least the variables:
%       - bepcode       An ID, unique for each lab test
%       - refwaarde     The variable that contains the one or double sided
%                       reference interval
%       - uitslag       The outcome of the lab test
%
% The data is standardized using the Z-score ( Z=(x-mean)/standard
% deviation)
%
% As output, the code generates adds a new variable to the LabResults field
% with the new standardized value.
%

%% Load datasets
load([subfold,'/',d,'/ImportedTables.mat'],'LabResults', ...
    'total_patientsID','GeneralData')

%% Create columns to be added to LabResults table

col_2_add = array2table(NaN(height(LabResults),3));
col_2_add.Properties.VariableNames = {'RefLow','RefUp','StandardizedResult'};

LabResults = [LabResults col_2_add];

clearvars col_2_add

%% Split the field 'refwaarde' to two seperate fields for lower and upper
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

%% Compute standardized scores
% Update A
A = LabResults(patientidx,:);
% Z-transformation
x = table2array(A(:,{'uitslagNUM' 'RefLow' 'RefUp'}));

z = NaN(length(x),1);

for i=1:length(x)
    if and(~isnan(x(i,2)),~isnan(x(i,3)))
        % Both values are available
        s = (x(i,3) - x(i,2))/4; % sigma = (RefUp - RefLow)/4
        m = x(i,2) + SD*s;        % mu = RefLow + 2s
        z(i,1) = (x(i,1)-m)/s;
    elseif and(~isnan(x(i,2)),isnan(x(i,3)))
        % RefUp not available, one sided/tail 
        z(i,1) = (x(i,1)/x(i,2))*4-6;
    end
end

Standardized = array2table(z);
Standardized.Properties.VariableNames = {'StandardizedResult'};

% Allocate in original dataset
LabResults(patientidx,{'StandardizedResult'}) = Standardized;
clearvars StandardizedResult Standardized s m z
end

clear j i

BKR00b = LabResults(strcmp(LabResults.bepcode,'BKR000'),{'PatientCode',...
    'uitslagNUM','RefLow','RefUp','days_from_ok'});
ASAT = LabResults(strcmp(LabResults.bepcode,'BAS002'),{'PatientCode',...
    'days_from_ok','uitslagNUM','RefLow','RefUp'});
ALAT = LabResults(strcmp(LabResults.bepcode,'BAL015'),{'PatientCode',...
    'days_from_ok','uitslagNUM','RefLow','RefUp'}); 
Albumin = LabResults(strcmp(LabResults.bepcode,'BAL002'),{'PatientCode',...
    'uitslagNUM','days_from_ok','RefLow','RefUp'});
Calcium = LabResults(strcmp(LabResults.bepcode,'BCA006'),{'PatientCode',...
    'uitslagNUM','days_from_ok','RefLow','RefUp'});

%% Standarize Age

ZAge = array2table(NaN(size(total_patientsID,1),2));
ZAge.Properties.VariableNames = {'PatientCode','Age'};
[D,PD] = allfitdist(GeneralData.Age);

mAge = NaN(3,1);
sAge = NaN(3,1);

DAge = D; PDAge = PD;
    
for i=1:3
    s = std(PD{i}); sAge(i,1) = s;
    m = mean(PD{i}); mAge(i,1) = m;
    figure('Name',['Age ',D(i).DistName])
    nbins = max(min(length(GeneralData.Age)./10,100),50);
    xi = linspace(min(GeneralData.Age),max(GeneralData.Age),nbins);
    dx = mean(diff(xi));
    xi2 = linspace(min(GeneralData.Age),max(GeneralData.Age),nbins*10)';
    fi = histc(GeneralData.Age,xi-dx);
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
    xlabel('Age','FontSize',18,'FontWeight','bold') % x-axis label
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

        print([subfold,'/',d,'/Plots/Age',D(i).DistName],'-deps')
        print([subfold,'/',d,'/Plots/Age',D(i).DistName,'c'],'-depsc')
%         print([subfold,'/',d,'/Plots/Age',D(i).DistName], '-dpdf','-bestfit')
end

close all
clearvars ax b1 D dx fi fig h1 h2 h3 h4 h5 i m nbins p1 PD s xi xi2 xL ys

% Standarize values
s = ((mAge(1,1)+2*sAge(1,1)) - (mAge(1,1)-2*sAge(1,1)))/4; % sigma = (RefUp - RefLow)/4
m = mAge(1,1);       % mu = RefLow + 2s
for pp=1:size(total_patientsID,1)
    ZAge(pp,{'PatientCode'}) = array2table(GeneralData.PatientNr(pp));
    ZAge(pp,{'Age'}) = array2table((table2array(GeneralData(pp,{'Age'}))-m)/s);
end
clearvars s m pp
