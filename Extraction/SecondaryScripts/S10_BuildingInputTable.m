%% Create input table

close all

% load data

load([subfold,'/',d,'/Classifiers',char(string(year(datetime(date)))),...
    '+-3months.mat'],'Classifiers2yr','Classifiers1yr')

load([subfold,'/',d,'/ImportedTables.mat'],'ScreeningData', ...
    'total_patientsID', 'GeneralData')

load([subfold,'/',d,'/LabResultsTimeSeries.mat'], 'LabLastPatient')

load([subfold,'/',d,'/PhMsTimeSeries.mat'], 'PhMsLastPatient')

load([subfold,'/',d,'/PPOSTimeSeries.mat'], 'PPOSMedianPatient',...
    'PPOSLastPatient')

load([subfold,'/',d,'/PreOperativeValues.mat'], 'PreOpLast', 'PreOpMax',...
    'PreOpMean', 'PreOpMedian', 'PreOpMin', 'bmi_max')

load([subfold,'/',d,'/LabZt.mat'], 'ZAge')

%% Define options
if create_folder == true
    % Define folder name
    root = 'Input_tables';
    if bin_class == true
        foldername = sprintf('%s_(%s)',root, id_twl);
    end
    if multi_class == true
        foldername = sprintf('%s_(%s)',root, 'multi');
    end
    if reg_model == true
        foldername = sprintf('%s_(%s)',root, 'reg');
    end
    comp_foldername = [subfold,'/',d,'/',foldername];
    if isfolder([subfold,'/',d,'/',foldername]) == false
        mkdir(comp_foldername)
    end
end

% Format type .fig MATLAB default
if keep_fig == true
    fttype = '.png';
end
%% Extract TWL patients and discretisize

if bin_class == true
    % After a year
    Output1yr = Classifiers1yr(:,1:2);
    Output1yr(isnan(Output1yr.TWL),:)= [];

    figure
    histogram(Output1yr.TWL,'NumBins',10,'FaceColor', [1 0.2 0.2],'BinWidth',0.1)
    xlim([0 1]);
    if keep_fig == true
        ylabel('Frequency','FontSize',16,'FontWeight','bold')
        legend('Total Weight Loss %','FontSize',12);
        title('Distribution of TWL% after 1yr','FontSize',18)
        ax = gca; ax.FontSize = 12;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of TWL% after 1yr',fttype]);
    end
    
    idx_class0_1yr = [Output1yr.TWL] >= thr_twl; % at least 20% TWL, suggested 35%
    total_sucess = sum(idx_class0_1yr);
    total_non_sucess = size(Output1yr,1)- total_sucess;
    Output1yr.TWL(idx_class0_1yr) = 1;
    Output1yr.TWL(~idx_class0_1yr) = 0;

    % After 2 yrs
    Output2yr = Classifiers2yr(:,1:2);
    Output2yr(isnan(Output2yr.TWL),:)= [];

    figure
    histogram(Output2yr.TWL,'NumBins',10,'FaceColor', [1 0.2 0.2],'BinWidth',0.1)
    xlim([0 1]);
    if keep_fig == true
        ylabel('Frequency','FontSize',16,'FontWeight','bold')
        legend('Total Weight Loss %','FontSize',12);
        title('Distribution of TWL% after 2 yrs','FontSize',18)
        ax = gca; ax.FontSize = 12;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of TWL% after 2 yrs',fttype]);
    end
    
    idx_sucess_2yr = [Output2yr.TWL] >= thr_twl;
    total_sucess2 = sum(idx_sucess_2yr);
    total_non_sucess2 = size(Output2yr,1)- total_sucess2;
    Output2yr.TWL(idx_sucess_2yr) = 1;
    Output2yr.TWL(~idx_sucess_2yr) = 0;
end

if multi_class == true
    % After a year
    Output1yr = Classifiers1yr(:,1:2);
    Output1yr(isnan(Output1yr.TWL),:)= [];

    figure
    title('Distribution of TWL% after 1yr')
    histogram(Output1yr.TWL,'NumBins',10,'FaceColor', [1 0.2 0.2],'BinWidth',0.1)
    xlim([0 1]);
    if keep_fig == true
        %savefig('Distribution of TWL% after 1yr.png')
        saveas(gcf,[comp_foldername '/Distribution of TWL% after 1yr',fttype]);
    end
    
    idx_class0_1yr = [Output1yr.TWL] <= class0;
    total_class1_1yr = sum(idx_class0_1yr);
    idx_class1_1yr = arrayfun(@(x) x > class0 & x <= class1, [Output1yr.TWL]);
    total_class2_1yr = sum(idx_class1_1yr);
    idx_class2_1yr = arrayfun(@(x) x > class1 & x <= class2, [Output1yr.TWL]);
    total_class3_1yr = sum(idx_class2_1yr);
    idx_class3_1yr = [Output1yr.TWL] > class2;
    total_class4_1yr = sum(idx_class3_1yr);
    Output1yr.TWL(idx_class0_1yr) = 1;
    Output1yr.TWL(idx_class1_1yr) = 2;
    Output1yr.TWL(idx_class2_1yr) = 3;
    Output1yr.TWL(idx_class3_1yr) = 4;
    
    Output1yr.TWL = double(ordinal([Output1yr.TWL],{'0','1','2','3'},[],[0,class0,class1,class2,Inf]));
    total_class1_1yr = sum([Output1yr.TWL]==1);
    total_class2_1yr = sum([Output1yr.TWL]==2);
    total_class3_1yr = sum([Output1yr.TWL]==3);
    total_class4_1yr = sum([Output1yr.TWL]==4);
    
    % After 2 yrs
    Output2yr = Classifiers2yr(:,1:2);
    Output2yr(isnan(Output2yr.TWL),:)= [];

    figure
    title('Distribution of TWL% after 2 yrs')
    histogram(Output2yr.TWL,'NumBins',10,'FaceColor', [1 0.2 0.2],'BinWidth',0.1)
    xlim([0 1]);
    if keep_fig == true
        saveas(gcf,[comp_foldername '/Distribution of TWL% after 2 yrs',fttype]);
    end
    
    Output2yr.TWL = double(ordinal([Output2yr.TWL],{'0','1','2','3'},[],[0,class0,class1,class2,Inf]));
    total_class1_2yr = sum([Output2yr.TWL]==1);
    total_class2_2yr = sum([Output2yr.TWL]==2);
    total_class3_2yr = sum([Output2yr.TWL]==3);
    total_class4_2yr = sum([Output2yr.TWL]==4);
    
    clearvars idx_class0_1yr idx_class0_2yr idx_class1_1yr idx_class1_2yr ...
        idx_class2_1yr idx_class2_2yr idx_class3_1yr idx_class3_2yr id_twl
end

if reg_model == true
    % After a year
    Output1yr = Classifiers1yr(:,1:2);
    Output1yr(isnan(Output1yr.TWL),:)= [];

    figure
    title('Distribution of TWL% after 1yr')
    histogram(Output1yr.TWL,'NumBins',10,'FaceColor', [1 0.2 0.2],'BinWidth',0.1)
    xlim([0 inf]);
    if keep_fig == true
        % savefig('Distribution of TWL% after 1yr.png')
        saveas(gcf,[comp_foldername '/Distribution of TWL% after 1yr',fttype]);
    end
    
    % After 2 yrs
    Output2yr = Classifiers2yr(:,1:2);
    Output2yr(isnan(Output2yr.TWL),:)= [];

    figure
    title('Distribution of TWL% after 2 yrs')
    histogram(Output2yr.TWL,'NumBins',10,'FaceColor', [1 0.2 0.2],'BinWidth',0.1)
    xlim([0 1]);
    if keep_fig == true
        saveas(gcf,[comp_foldername '/Distribution of TWL% after 2 yrs',fttype]);
    end
end

%% Add Input data
Lab = size(LabLastPatient,2);
Lab_names = LabLastPatient.Properties.VariableNames;

PhMs_names = PhMsLastPatient.Properties.VariableNames(end-1:end);
PhMs = length(PhMs_names);

PPOS_names = [];
PPOS = 0;

GData = size(GeneralData,2)-2;
GData_names = GeneralData.Properties.VariableNames(3:end);

Screen_names = {'roken','rokenpy','diabet','hypert','dyslip'};
Screen = length(Screen_names);

PreOp = 3;
PreOp_names = {'BMI','WC','TBFM'};

input_size_columns = Lab+PhMs+PPOS+GData+Screen+PreOp;

var_names = horzcat(Lab_names,PhMs_names,PPOS_names,GData_names,...
    Screen_names,PreOp_names);

clearvars Lab PhMs PPOS SuType GData Screen PreOp Lab_names PhMs_names ...
    PPOS_names GData_names PreOp_names;

%% Pre-allocate
input_1yr = array2table(NaN(size(Output1yr,1),input_size_columns+6));
input_1yr.Properties.VariableNames = horzcat(var_names,...
    Classifiers1yr.Properties.VariableNames{2:5},{'Age_','Metabolic'});

input_2yr = array2table(NaN(size(Output2yr,1),input_size_columns+6));
input_2yr.Properties.VariableNames = horzcat(var_names,...
    Classifiers2yr.Properties.VariableNames{2:5},{'Age_','Metabolic'});

%% Allocate data after 1 year
p = size(var_names,2)-1;

for i=1:size(Output1yr,1)
    patientID = Output1yr.PatientCode(i);
    input_1yr(i,{'PatientCode'})= array2table(patientID);
    input_1yr(i,{'TWL'}) = array2table(Output1yr.TWL(i));
    input_1yr(i,{'WL'}) = array2table(Classifiers1yr.WL(i));
    input_1yr(i,{'EBMIL'}) = array2table(Classifiers1yr.EBMIL(i));
    input_1yr(i,{'EBWL'}) = array2table(Classifiers1yr.EBWL(i));
    
    % Find indexes in other tables
    idx_lab = [LabLastPatient.PatientCode] == patientID;
    idx_PhMs = (cell2mat(PhMsLastPatient.PatientCode)) == patientID;
    idx_PPOS = [PPOSLastPatient.PatientCode] == patientID;
    idx_SBP = (cell2mat(PhMsLastPatient.PatientCode)) == patientID;
    idx_DBP = (cell2mat(PhMsLastPatient.PatientCode)) == patientID;
    idx_GData = [GeneralData.PatientNr] == patientID;
    idx_Screen = [ScreeningData.PatientCode] == patientID;
    idx_PreOp = [PreOpLast.PatientCode] == patientID;
    idx_bmi_max = [bmi_max.PatientCode] == patientID;
    
    input_1yr(i,{'Metabolic'}) = bmi_max(idx_bmi_max,{'Max'});
    
    for lab=1:size(LabLastPatient(:,2:end),2)
        lab_name = LabLastPatient.Properties.VariableNames{lab+1};
        input_1yr(i,lab_name) = LabLastPatient(idx_lab,lab_name);
    end
    clearvars lab lab_name

    input_1yr(i,{'Systolic'}) = array2table(cell2mat(PhMsLastPatient.Systolic(idx_SBP)));
    
    mintension_patient = cell2mat(PhMsLastPatient.Diastolic(idx_DBP));
    if mintension_patient == 0
       input_1yr(i,{'Diastolic'}) = array2table(NaN);
    else
        input_1yr(i,{'Diastolic'}) = array2table(mintension_patient);
    end
    
    input_1yr(i,'Age') = ZAge(idx_GData,{'Age'});
    input_1yr(i,'Age_') = GeneralData(idx_GData,{'Age'});
    input_1yr(i,'Sex') = GeneralData(idx_GData,{'Sex'});
    
    for screen = 1:length(Screen_names)
        screen_name = Screen_names{screen};
        input_1yr(i,screen_name) = ScreeningData(idx_Screen,screen_name);
    end
    clearvars screen_name

    input_1yr(i,{'BMI'}) = PreOpLast(idx_PreOp,{'BMI'});
    input_1yr(i,{'WC'}) = PreOpLast(idx_PreOp,{'WaistCircum'});
    input_1yr(i,{'TBFM'}) = PreOpLast(idx_PreOp,{'TBFM'});
end

clearvars i

mat_1yr = table2array(input_1yr);

dsize = size(mat_1yr,1);

figure('Name','Missingness input 1yr')
imagesc(isnan(mat_1yr(:,2:end)));
if keep_fig == true
    title(['After gathering all selected variables n=', num2str(dsize)])
    xlabel('Variables','FontSize',16)
    ylabel('Yellow = missing, Blue = no missing','FontSize',16)
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/After gathering all selected variables for ', num2str(dsize),' patients',fttype]);
end
close all
%% Distribution of success class among gender, type of surgery, BMI and comorbidities (1yr)
if hist_dist == true
    conta = 0;
    % Gender
    male = input_1yr.Sex == 1;
    female = input_1yr.Sex == 0;
    
    figure('Name','Gender 1yr')
    histogram(input_1yr{:,'TWL'}(female),'BinMethod', 'integer',...
        'FaceColor', '#CF0047','EdgeColor','none','FaceAlpha',0.4);
        hold on
    histogram(input_1yr{:,'TWL'}(male),'BinMethod', 'integer',...
        'FaceColor', '#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    legend('female','male','FontSize',12)
    xlabel('Output variable (%TWL)','FontSize',16,...
        'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of TWL among gender (1yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among gender (1yr, n=', num2str(dsize),',',id_twl,'B)',fttype]);
    end

    figure('Name','Age by sex')
        hold on
    histogram(input_1yr{:,'Age_'}(female),'BinMethod', 'integer',...
        'FaceColor','#CF0047','FaceAlpha',0.4,'EdgeColor','none');
    histogram(input_1yr{:,'Age_'}(male),'BinMethod', 'integer',...
        'FaceColor','#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    hold off
    legend('female','male','FontSize',12)
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
    if keep_fig == true
        title(['Distribution of ages by sex (1yr, n=',...
        num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of ages (1 yr n=', num2str(dsize),',',id_twl,')',fttype]);
    end
    
    figure('Name','Normplot Age')
    h1 = normplot(input_1yr{:,'Age_'});
    h1(1).LineWidth = 2;h1(2).LineWidth = 2;h1(3).LineWidth = 2;
    h1(1).MarkerEdgeColor = '#E5CC00';
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age (1 year, n=', ...
        num2str(dsize),',',id_twl,')',fttype]);
    
    figure('Name','Normplots Age by sex')
    hold on
    hm = normplot(input_1yr{:,'Age_'}(male));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(input_1yr{:,'Age_'}(female));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Male','','','Female'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age by sex (1 yr, n=', ...
        num2str(dsize),',',id_twl,')',fttype]);
    
    % Null hypothesis: The mean strengths for the two populations are equal.
    % A p-value less than the significance level indicates that you 
    % can reject the null hypothesis. In other words, the sample provides 
    % sufficient evidence to conclude that the population means are 
    % different. Below is the output for the analysis.
    [httest2,pttest2,cittest2,statsttest2] = ttest2((input_1yr{:,'Age_'}(male)),...
    (input_1yr{:,'Age_'}(female)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((input_1yr{:,'Age_'}(male)),...
    (input_1yr{:,'Age_'}(female)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.A.first.tstudent(conta+1).fieldName = 'Age by sex';
    stat.A.first.tstudent(conta+1).h_ttest = httest2;
    stat.A.first.tstudent(conta+1).pval_ttest = pttest2;
    stat.A.first.tstudent(conta+1).ci_ttest = cittest2;
    stat.A.first.tstudent(conta+1).tstat = statsttest2.tstat;
    stat.A.first.tstudent(conta+1).df = statsttest2.df;
    stat.A.first.tstudent(conta+1).sd = statsttest2.sd;
    stat.A.first.tstudent(conta+1).tstudent = veredict1;
    stat.A.first.tstudent(conta+1).male = male;
    stat.A.first.tstudent(conta+1).female = female;
    stat.A.first.tstudent(conta+1).data = (input_1yr(:,{'PatientCode','Sex','Age_'}));
    
    
    stat.A.first.Ftest(conta+1).fieldName = 'Age by sex';
    stat.A.first.Ftest(conta+1).fstat = stats_var.fstat;
    stat.A.first.Ftest(conta+1).df1 = stats_var.df1;
    stat.A.first.Ftest(conta+1).df2 = stats_var.df2;
    stat.A.first.Ftest(conta+1).h_ftest = h_var;
    stat.A.first.Ftest(conta+1).pval_ftest = p_var;
    stat.A.first.Ftest(conta+1).ci_ftest = ci_var;
    stat.A.first.Ftest(conta+1).ftest = veredict2;
    stat.A.first.Ftest(conta+1).male = male;
    stat.A.first.Ftest(conta+1).female = female;
    stat.A.first.Ftest(conta+1).data = (input_1yr(:,{'PatientCode','Sex','Age_'}));
    
    tbl_cross = array2table(crosstab(female,(input_1yr{:,'TWL'}))); % 1 = female
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the 2 distinct GENDERS.
    tbl_cross.Properties.RowNames = {'male','female'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association gender & TWL ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association gender & TWL';
    end
    
    stat.A.first.Fisher(conta+1).fieldName = 'Gender by class';
    stat.A.first.Fisher(conta+1).h = h;
    stat.A.first.Fisher(conta+1).pval = p;
    stat.A.first.Fisher(conta+1).OddsRatio = stats.OddsRatio;
    stat.A.first.Fisher(conta+1).ci = stats.ConfidenceInterval;
    stat.A.first.Fisher(conta+1).Fisher = veredict3;
    stat.A.first.Fisher(conta+1).contingencyTable = tbl_cross;
    stat.A.first.Fisher(conta+1).yes = 'TWL';
    stat.A.first.Fisher(conta+1).no = 'TWL';
    stat.A.first.Fisher(conta+1).data = (input_1yr(:,{'PatientCode','TWL','Age_'}));
    
    clearvars male female
    
    % BMI
    success = input_1yr.TWL == 1;
    fail = input_1yr.TWL == 0;
    
    figure('Name','BMI 1st yr')
    histogram(table2array(input_1yr(success,{'BMI'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(input_1yr(fail,{'BMI'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('BMI (kg/m^2)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class among BMI (1yr, n=', num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among BMI min (1yr, n=', num2str(dsize),',',id_twl,')',fttype]);
    end
        
    figure('Name','Normplots BMI')
    hold on
    hm = normplot(input_1yr{:,'BMI'}(success));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(input_1yr{:,'BMI'}(fail));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Success','','','Fail'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of BMI by class (1 yr, n=', ...
        num2str(dsize),',',id_twl,')',fttype]);
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((input_1yr{:,'BMI'}(success)),...
    (input_1yr{:,'BMI'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((input_1yr{:,'BMI'}(success)),...
    (input_1yr{:,'BMI'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.A.first.tstudent(conta+2).fieldName = 'BMI';
    stat.A.first.tstudent(conta+2).h_ttest = httest2;
    stat.A.first.tstudent(conta+2).pval_ttest = pttest2;
    stat.A.first.tstudent(conta+2).ci_ttest = cittest2;
    stat.A.first.tstudent(conta+2).tstat = statsttest2.tstat;
    stat.A.first.tstudent(conta+2).df = statsttest2.df;
    stat.A.first.tstudent(conta+2).sd = statsttest2.sd;
    stat.A.first.tstudent(conta+2).tstudent = veredict1;
    stat.A.first.tstudent(conta+2).sucess = success;
    stat.A.first.tstudent(conta+2).fail = fail;
    stat.A.first.tstudent(conta+2).data = (input_1yr(:,{'PatientCode','TWL','BMI'}));
    
    [pWhitney,hWhitney,statsWhitney] = ranksum((input_1yr{:,'BMI'}(success)),...
    (input_1yr{:,'BMI'}(fail)));
    if (hWhitney == 1)
        % rejects the null hypothesis
        veredictW = 'Diff sig (median)'; 
    else
        veredictW = 'No diff sig (median)';
    end
    
    stat.A.first.Whitney(conta+1).fieldName = 'BMI';
    stat.A.first.Whitney(conta+1).h = hWhitney;
    stat.A.first.Whitney(conta+1).p = pWhitney;
    stat.A.first.Whitney(conta+1).ranksum = statsWhitney.ranksum;
    stat.A.first.Whitney(conta+1).Whitney = veredictW;
    stat.A.first.Whitney(conta+1).sucess = success;
    stat.A.first.Whitney(conta+1).fail = fail;
    stat.A.first.Whitney(conta+1).data = (input_1yr(:,{'PatientCode','TWL','BMI'}));
    
    stat.A.first.Ftest(conta+2).fieldName = 'BMI';
    stat.A.first.Ftest(conta+2).fstat = stats_var.fstat;
    stat.A.first.Ftest(conta+2).df1 = stats_var.df1;
    stat.A.first.Ftest(conta+2).df2 = stats_var.df2;
    stat.A.first.Ftest(conta+2).h_ftest = h_var;
    stat.A.first.Ftest(conta+2).pval_ftest = p_var;
    stat.A.first.Ftest(conta+2).ci_ftest = ci_var;
    stat.A.first.Ftest(conta+2).ftest = veredict2;
    stat.A.first.Ftest(conta+2).sucess = success;
    stat.A.first.Ftest(conta+2).fail = fail;
    stat.A.first.Ftest(conta+2).data = (input_1yr(:,{'PatientCode','TWL','BMI'}));
    
    % Age
    figure('Name','Age 1 yr')
    histogram(table2array(input_1yr(success,{'Age_'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(input_1yr(fail,{'Age_'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
        hold off
        if keep_fig == true
            title(['Distribution of ages among classes (1yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of ages among classes (1yr, n=', num2str(dsize),',',id_twl,')',fttype]);
        end
    
    figure('Name','Normplot Age (class)')
    hold on
    hm = normplot(input_1yr{:,'Age_'}(success));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(input_1yr{:,'Age_'}(fail));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Success','','','Fail'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of TWL by sex (1 yr, n=', ...
        num2str(dsize),',',id_twl,')',fttype]);
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((input_1yr{:,'Age_'}(success)),...
    (input_1yr{:,'Age_'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((input_1yr{:,'Age_'}(success)),...
    (input_1yr{:,'Age_'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.A.first.tstudent(conta+3).fieldName = 'Age and Class';
    stat.A.first.tstudent(conta+3).h_ttest = httest2;
    stat.A.first.tstudent(conta+3).pval_ttest = pttest2;
    stat.A.first.tstudent(conta+3).ci_ttest = cittest2;
    stat.A.first.tstudent(conta+3).tstat = statsttest2.tstat;
    stat.A.first.tstudent(conta+3).df = statsttest2.df;
    stat.A.first.tstudent(conta+3).sd = statsttest2.sd;
    stat.A.first.tstudent(conta+3).tstudent = veredict1;
    stat.A.first.tstudent(conta+3).sucess = success;
    stat.A.first.tstudent(conta+3).fail = fail;
    stat.A.first.tstudent(conta+3).data = (input_1yr(:,{'PatientCode','TWL','BMI'}));
    
    stat.A.first.Ftest(conta+3).fieldName = 'BMI';
    stat.A.first.Ftest(conta+3).fstat = stats_var.fstat;
    stat.A.first.Ftest(conta+3).df1 = stats_var.df1;
    stat.A.first.Ftest(conta+3).df2 = stats_var.df2;
    stat.A.first.Ftest(conta+3).h_ftest = h_var;
    stat.A.first.Ftest(conta+3).pval_ftest = p_var;
    stat.A.first.Ftest(conta+3).ci_ftest = ci_var;
    stat.A.first.Ftest(conta+3).ftest = veredict2;
    stat.A.first.Ftest(conta+3).sucess = success;
    stat.A.first.Ftest(conta+3).fail = fail;
    stat.A.first.Ftest(conta+3).data = (input_1yr(:,{'PatientCode','TWL','Age_'}));
    
    clearvars success fail
    
    % Dyslipidemia
    yes = input_1yr.dyslip == 1;
    no = input_1yr.dyslip == 0;
    
    figure('Name','Dyslipidemia 1 yr')
    histogram(table2array(input_1yr(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(input_1yr(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,...
        'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and dyslipedemia (1yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and dyslipedemia (1yr, n=', num2str(dsize),',',id_twl,')',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(input_1yr{:,'TWL'}))); % 1 = dyslipidemia
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_dyslip','dyslip'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association dyslipidemia&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association dyslipidemia&TWL';
    end
    
    stat.A.first.Fisher(conta+2).fieldName = 'Dyslipidemia by class';
    stat.A.first.Fisher(conta+2).h = h;
    stat.A.first.Fisher(conta+2).pval = p;
    stat.A.first.Fisher(conta+2).OddsRatio = stats.OddsRatio;
    stat.A.first.Fisher(conta+2).ci = stats.ConfidenceInterval;
    stat.A.first.Fisher(conta+2).Fisher = veredict3;
    stat.A.first.Fisher(conta+2).yes = yes;
    stat.A.first.Fisher(conta+2).no = no;
    stat.A.first.Fisher(conta+2).contingencyTable = tbl_cross;
    stat.A.first.Fisher(conta+2).data = (input_1yr(:,{'PatientCode','TWL','dyslip'}));
    
    clearvars yes no
    
    % Diabetes
    yes = input_1yr.diabet == 1;
    no = input_1yr.diabet == 0;
    
    figure('Name','Diabetes 1 yr')
    histogram(table2array(input_1yr(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(input_1yr(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
     if keep_fig == true
        title(['Distribution of success class and diabetes (1yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and diabetes (1yr, n=', num2str(dsize),',',id_twl,')',fttype]);
     end
     
     tbl_cross = array2table(crosstab(yes,(input_1yr{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_diabetes','diabetes'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association diabetes&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association diabetes&TWL';
    end
    
    stat.A.first.Fisher(conta+3).fieldName = 'Diabetes by class';
    stat.A.first.Fisher(conta+3).h = h;
    stat.A.first.Fisher(conta+3).pval = p;
    stat.A.first.Fisher(conta+3).OddsRatio = stats.OddsRatio;
    stat.A.first.Fisher(conta+3).ci = stats.ConfidenceInterval;
    stat.A.first.Fisher(conta+3).Fisher = veredict3;
    stat.A.first.Fisher(conta+3).yes = yes;
    stat.A.first.Fisher(conta+3).no = no;
    stat.A.first.Fisher(conta+3).contingencyTable = tbl_cross;
    stat.A.first.Fisher(conta+3).data = (input_1yr(:,{'PatientCode','TWL','diabet'}));
    
    clearvars yes no
    
    % Hypertension
    yes = input_1yr.hypert == 1;
    no = input_1yr.hypert == 0;
    
    figure('Name','Hypertension 1 yrs')
    histogram(table2array(input_1yr(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(input_1yr(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and hypertension (1yr, n=',num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and hypertension (1yr, n=',num2str(dsize),',',id_twl,')',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(input_1yr{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_hypert','hypert'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association hypertension&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association hypertension&TWL';
    end
    
    stat.A.first.Fisher(conta+4).fieldName = 'Hypertension by class';
    stat.A.first.Fisher(conta+4).h = h;
    stat.A.first.Fisher(conta+4).pval = p;
    stat.A.first.Fisher(conta+4).OddsRatio = stats.OddsRatio;
    stat.A.first.Fisher(conta+4).ci = stats.ConfidenceInterval;
    stat.A.first.Fisher(conta+4).Fisher = veredict3;
    stat.A.first.Fisher(conta+4).yes = yes;
    stat.A.first.Fisher(conta+4).no = no;
    stat.A.first.Fisher(conta+4).contingencyTable = tbl_cross;
    stat.A.first.Fisher(conta+4).data = (input_1yr(:,{'PatientCode','TWL','hypert'}));
    
    clearvars yes no
end

close all
%% Allocate data after 2 year

for i=1:size(Output2yr,1)
    patientID = Output2yr.PatientCode(i);
    input_2yr(i,{'PatientCode'})= array2table(patientID);
    input_2yr(i,{'TWL'}) = array2table(Output2yr.TWL(i));
    input_2yr(i,{'WL'}) = array2table(Classifiers2yr.WL(i));
    input_2yr(i,{'EBMIL'}) = array2table(Classifiers2yr.EBMIL(i));
    input_2yr(i,{'EBWL'}) = array2table(Classifiers2yr.EBWL(i));
    
    % Find indexes in other tables
    idx_lab = [LabLastPatient.PatientCode] == patientID;
    idx_PhMs = (cell2mat(PhMsLastPatient.PatientCode)) == patientID;
    idx_PPOS = [PPOSLastPatient.PatientCode] == patientID;
    idx_SBP = (cell2mat(PhMsLastPatient.PatientCode)) == patientID;
    idx_DBP = (cell2mat(PhMsLastPatient.PatientCode)) == patientID;
    idx_GData = [GeneralData.PatientNr] == patientID;
    idx_Screen = [ScreeningData.PatientCode] == patientID;
    idx_PreOp = [PreOpLast.PatientCode] == patientID;
    idx_bmi_max = [bmi_max.PatientCode] == patientID;
    
    input_2yr(i,{'Metabolic'}) = bmi_max(idx_bmi_max,{'Max'});
    
    for lab=1:size(LabLastPatient(:,2:end),2)
        lab_name = LabLastPatient.Properties.VariableNames{lab+1};
        input_2yr(i,lab_name) = LabLastPatient(idx_lab,lab_name);
    end
    clearvars lab lab_name
    
    input_2yr(i,{'Systolic'}) = array2table(cell2mat(PhMsLastPatient.Systolic(idx_SBP)));
    
    mintension_patient = cell2mat(PhMsLastPatient.Diastolic(idx_DBP));
    if mintension_patient == 0
       input_2yr(i,{'Diastolic'}) = array2table(NaN);
    else
        input_2yr(i,{'Diastolic'}) = array2table(mintension_patient);
    end
    clearvars mintension_patient
    
    input_2yr(i,'Age') = ZAge(idx_GData,{'Age'});
    input_2yr(i,'Sex') = GeneralData(idx_GData,{'Sex'});
    input_2yr(i,'Age_') = GeneralData(idx_GData,{'Age'});

    for screen = 1:length(Screen_names)
        screen_name = Screen_names{screen};
        tmp = ScreeningData(idx_Screen,screen_name);
        if isempty(tmp) == 1
            input_2yr(i,screen_name) = array2table(NaN);
        else
            input_2yr(i,screen_name) = tmp;
        end
    end
    
    clearvars screen_name tmp
    
    input_2yr(i,{'BMI'}) = PreOpLast(idx_PreOp,{'BMI'});
    input_2yr(i,{'WC'}) = PreOpLast(idx_PreOp,{'WaistCircum'});
    input_2yr(i,{'TBFM'}) = PreOpLast(idx_PreOp,{'TBFM'});
end

clearvars i
    
mat_2yr = table2array(input_2yr);

dsize = size(mat_2yr,1);

figure('Name','Missingness input 2 yrs')
imagesc(isnan(mat_2yr(:,2:end)));
if keep_fig == true
    title(['After gathering all selected variables n=', num2str(dsize)])
    ylabel('Yellow = missing, Blue = no missing','FontSize',16)
    xlabel('Variables','FontSize',16)
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/After gathering all selected variables for ', num2str(dsize) ,' patients',fttype]);
end

close all

%% Distribution of success class among gender, type of surgery, BMI and comorbidities (2yrs)

if hist_dist == true
     if bin_class == true
        % Gender
        male = input_2yr.Sex == 1;
        female = input_2yr.Sex == 0;

        figure('Name','Gender 2 yrs')
        histogram(input_2yr{:,'TWL'}(female),'BinMethod', 'integer', ...
            'FaceColor', '#CF0047','EdgeColor','none','FaceAlpha',0.4);
            hold on
        histogram(input_2yr{:,'TWL'}(male),'BinMethod', 'integer',...
            'FaceColor', '#0097CD','FaceAlpha',0.4,'EdgeColor','none');
        legend('female','male','FontSize',12)
        xlabel('Output variable (%TWL)','FontSize',16,...
            'FontWeight','bold')
            hold off
        if keep_fig == true
            title(['Distribution of success class among gender (2yr, n=',...
                num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax= gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 25; %Width
            fig.Position(4) = 18.7; %Height
            saveas(gcf,[comp_foldername '/Distribution of success class among gender (2yr, n=', num2str(dsize),',',id_twl,')',fttype]);
        end
        
        figure('Name','Age by sex')
            hold on
        histogram(input_2yr{:,'Age_'}(female),'BinMethod', 'integer',...
            'FaceColor','#CF0047','FaceAlpha',0.4,'EdgeColor','none');
        histogram(input_2yr{:,'Age_'}(male),'BinMethod', 'integer',...
            'FaceColor','#0097CD','FaceAlpha',0.4,'EdgeColor','none');
        hold off
        legend('female','male','FontSize',12)
        xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
        if keep_fig == true
            title(['Distribution of ages by sex (2yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of ages (2 yr n=', num2str(dsize),',',id_twl,'B)',fttype]);
        end
        
        figure('Name','Normplot Age')
        h1 = normplot(input_2yr{:,'Age_'});
        h1(1).LineWidth = 2;h1(2).LineWidth = 2;h1(3).LineWidth = 2;
        h1(1).MarkerEdgeColor = '#E5CC00';
        xlabel('Age at surgery (2 year)','FontSize',16,'FontWeight','bold')
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Normplot of Age (2 year, n=', ...
            num2str(dsize),',',id_twl,')',fttype]);

        figure('Name','Normplots Age by sex')
        hold on
        hm = normplot(input_2yr{:,'Age_'}(male));
        hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
        hm(1).MarkerEdgeColor = '#0097CD';
        hf = normplot(input_2yr{:,'Age_'}(female));
        hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
        hf(1).MarkerEdgeColor = '#CF0047';
        hold off
        legend({'','','Male','','','Female'},'Location','southeast')
        xlabel('Age at surgery (2 year)','FontSize',16,'FontWeight','bold')
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Normplot of Age by sex (2 yr, n=', ...
            num2str(dsize),',',id_twl,')',fttype]);
        
        % Null hypothesis: The mean strengths for the two populations are equal.
        % A p-value less than the significance level indicates that you 
        % can reject the null hypothesis. In other words, the sample provides 
        % sufficient evidence to conclude that the population means are 
        % different. Below is the output for the analysis.
        [httest2,pttest2,cittest2,statsttest2] = ttest2((input_2yr{:,'Age_'}(male)),...
        (input_2yr{:,'Age_'}(female)));
        if (httest2 == 1)
            % rejects the null hypothesis
            veredict1 = 'Diff sig (mean)'; 
        else
            veredict1 = 'No diff sig (means)';
        end
        [h_var,p_var,ci_var,stats_var] = vartest2((input_2yr{:,'Age_'}(male)),...
        (input_2yr{:,'Age_'}(female)));
        if (h_var == 1)
            % rejects the null hypothesis
            veredict2 = 'Diff sig (var)'; 
        else
            veredict2 = 'No diff sig (var)';
        end

        stat.A.second.tstudent(conta+1).fieldName = 'Age by sex';
        stat.A.second.tstudent(conta+1).h_ttest = httest2;
        stat.A.second.tstudent(conta+1).pval_ttest = pttest2;
        stat.A.second.tstudent(conta+1).ci_ttest = cittest2;
        stat.A.second.tstudent(conta+1).tstat = statsttest2.tstat;
        stat.A.second.tstudent(conta+1).df = statsttest2.df;
        stat.A.second.tstudent(conta+1).sd = statsttest2.sd;
        stat.A.second.tstudent(conta+1).tstudent = veredict1;
        stat.A.second.tstudent(conta+1).male = male;
        stat.A.second.tstudent(conta+1).female = female;
        stat.A.second.tstudent(conta+1).data = (input_2yr(:,{'PatientCode','Sex','Age_'}));


        stat.A.second.Ftest(conta+1).fieldName = 'Age by sex';
        stat.A.second.Ftest(conta+1).fstat = stats_var.fstat;
        stat.A.second.Ftest(conta+1).df1 = stats_var.df1;
        stat.A.second.Ftest(conta+1).df2 = stats_var.df2;
        stat.A.second.Ftest(conta+1).h_ftest = h_var;
        stat.A.second.Ftest(conta+1).pval_ftest = p_var;
        stat.A.second.Ftest(conta+1).ci_ftest = ci_var;
        stat.A.second.Ftest(conta+1).ftest = veredict2;
        stat.A.second.Ftest(conta+1).male = male;
        stat.A.second.Ftest(conta+1).female = female;
        stat.A.second.Ftest(conta+1).data = (input_2yr(:,{'PatientCode','Sex','Age_'}));

        tbl_cross = array2table(crosstab(female,(input_2yr{:,'TWL'}))); % 1 = female
        %The rows in table correspond to the 2 distinct values in TWL, and
        %the columns correspond to the 2 distinct GENDERS.
        tbl_cross.Properties.RowNames = {'male','female'}; 
        tbl_cross.Properties.VariableNames = {'fail','success'};

        [h,p,stats] = fishertest(tbl_cross);
        if (h == 1)
            % Rejects null hypothesis
            veredict3 = ['Association gender & TWL ',...
                tbl_cross.Properties.RowNames{2}, ' have about ',...
                num2str(stats.OddsRatio),' greater odds of being ', ...
                tbl_cross.Properties.VariableNames{2}];
        else
            veredict3 = 'No association gender & TWL';
        end

        stat.A.second.Fisher(conta+1).fieldName = 'Gender by class';
        stat.A.second.Fisher(conta+1).h = h;
        stat.A.second.Fisher(conta+1).pval = p;
        stat.A.second.Fisher(conta+1).OddsRatio = stats.OddsRatio;
        stat.A.second.Fisher(conta+1).ci = stats.ConfidenceInterval;
        stat.A.second.Fisher(conta+1).Fisher = veredict3;
        stat.A.second.Fisher(conta+1).contingencyTable = tbl_cross;
        stat.A.second.Fisher(conta+1).yes = 'TWL';
        stat.A.second.Fisher(conta+1).no = 'TWL';
        stat.A.second.Fisher(conta+1).data = (input_2yr(:,{'PatientCode','TWL','Age_'}));
    
        clearvars male female

        % BMI
        success = input_2yr.TWL == 1;
        fail = input_2yr.TWL == 0;
        
        figure('Name','BMI 2nd year')
        histogram(table2array(input_2yr(success,{'BMI'})),...
            'BinMethod', 'integer', 'FaceColor', '#E19B00',...
            'EdgeColor','none')
            hold on
        histogram(table2array(input_2yr(fail,{'BMI'})),...
            'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
            'EdgeColor','none');
        legend('success','failure')
        xlabel('BMI (kg/m^2)','FontSize',16,'FontWeight','bold')
            hold off
        if keep_fig == true
            title(['Distribution of success class among BMI (2yrs, n=',...
                num2str(dsize),',',id_twl,')'])
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 25; %Width
            fig.Position(4) = 18.7; %Height
            saveas(gcf,[comp_foldername '/Distribution of success class among BMI min (2yrs, n=',num2str(dsize),',',id_twl,')',fttype]);
        end
        
        figure('Name','Normplots BMI')
        hold on
        hm = normplot(input_2yr{:,'BMI'}(success));
        hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
        hm(1).MarkerEdgeColor = '#0097CD';
        hf = normplot(input_2yr{:,'BMI'}(fail));
        hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
        hf(1).MarkerEdgeColor = '#CF0047';
        hold off
        legend({'','','Success','','','Fail'},'Location','southeast')
        xlabel('Age at surgery (2 year)','FontSize',16,'FontWeight','bold')
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Normplot of BMI by class (2 yr, n=', ...
            num2str(dsize),',',id_twl,')',fttype]);

        [httest2,pttest2,cittest2,statsttest2] = ttest2((input_2yr{:,'BMI'}(success)),...
        (input_2yr{:,'BMI'}(fail)));
        if (httest2 == 1)
            % rejects the null hypothesis
            veredict1 = 'Diff sig (mean)'; 
        else
            veredict1 = 'No diff sig (means)';
        end

        [h_var,p_var,ci_var,stats_var] = vartest2((input_2yr{:,'BMI'}(success)),...
        (input_2yr{:,'BMI'}(fail)));
        if (h_var == 1)
            % rejects the null hypothesis
            veredict2 = 'Diff sig (var)'; 
        else
            veredict2 = 'No diff sig (var)';
        end

        stat.A.second.tstudent(conta+2).fieldName = 'BMI';
        stat.A.second.tstudent(conta+2).h_ttest = httest2;
        stat.A.second.tstudent(conta+2).pval_ttest = pttest2;
        stat.A.second.tstudent(conta+2).ci_ttest = cittest2;
        stat.A.second.tstudent(conta+2).tstat = statsttest2.tstat;
        stat.A.second.tstudent(conta+2).df = statsttest2.df;
        stat.A.second.tstudent(conta+2).sd = statsttest2.sd;
        stat.A.second.tstudent(conta+2).tstudent = veredict1;
        stat.A.second.tstudent(conta+2).sucess = success;
        stat.A.second.tstudent(conta+2).fail = fail;
        stat.A.second.tstudent(conta+2).data = (input_2yr(:,{'PatientCode','TWL','BMI'}));
        
        [pWhitney,hWhitney,statsWhitney] = ranksum((input_2yr{:,'BMI'}(success)),...
        (input_2yr{:,'BMI'}(fail)));
        if (hWhitney == 1)
            % rejects the null hypothesis
            veredictW = 'Diff sig (median)';
        else
            veredictW = 'No diff sig (median)';
        end
        
        stat.A.second.Whitney(conta+1).fieldName = 'BMI';
        stat.A.second.Whitney(conta+1).h = hWhitney;
        stat.A.second.Whitney(conta+1).p = pWhitney;
        stat.A.second.Whitney(conta+1).ranksum = statsWhitney.ranksum;
        stat.A.second.Whitney(conta+1).Whitney = veredictW;
        stat.A.second.Whitney(conta+1).sucess = success;
        stat.A.second.Whitney(conta+1).fail = fail;
        stat.A.second.Whitney(conta+1).data = (input_1yr(:,{'PatientCode','TWL','BMI'}));
        
        stat.A.second.Ftest(conta+2).fieldName = 'BMI';
        stat.A.second.Ftest(conta+2).fstat = stats_var.fstat;
        stat.A.second.Ftest(conta+2).df1 = stats_var.df1;
        stat.A.second.Ftest(conta+2).df2 = stats_var.df2;
        stat.A.second.Ftest(conta+2).h_ftest = h_var;
        stat.A.second.Ftest(conta+2).pval_ftest = p_var;
        stat.A.second.Ftest(conta+2).ci_ftest = ci_var;
        stat.A.second.Ftest(conta+2).ftest = veredict2;
        stat.A.second.Ftest(conta+2).sucess = success;
        stat.A.second.Ftest(conta+2).fail = fail;
        stat.A.second.Ftest(conta+2).data = (input_2yr(:,{'PatientCode','TWL','BMI'}));

         % Age
        figure('Name','Age 2 yrs')
        histogram(table2array(input_2yr(success,{'Age_'})),...
            'BinMethod', 'integer', 'FaceColor', '#E19B00',...
            'EdgeColor','none')
        hold on
        histogram(table2array(input_2yr(fail,{'Age_'})),...
            'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
            'EdgeColor','none');
        legend('success','failure')
        xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
        hold off
        if keep_fig == true
            title(['Distribution of ages among classes (2yrs, n=',...
                num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of ages among classes (2yrs, n=',num2str(dsize),',',id_twl,')',fttype]);
        end
        
        [httest2,pttest2,cittest2,statsttest2] = ttest2((input_2yr{:,'Age_'}(success)),...
            (input_2yr{:,'Age_'}(fail)));
        if (httest2 == 1)
            % rejects the null hypothesis
            veredict1 = 'Diff sig (mean)';
        else
            veredict1 = 'No diff sig (means)';
        end
        
        [h_var,p_var,ci_var,stats_var] = vartest2((input_2yr{:,'Age_'}(success)),...
            (input_2yr{:,'Age_'}(fail)));
        if (h_var == 1)
            % rejects the null hypothesis
            veredict2 = 'Diff sig (var)';
        else
            veredict2 = 'No diff sig (var)';
        end
        
        stat.A.second.tstudent(conta+3).fieldName = 'Age and Class';
        stat.A.second.tstudent(conta+3).h_ttest = httest2;
        stat.A.second.tstudent(conta+3).pval_ttest = pttest2;
        stat.A.second.tstudent(conta+3).ci_ttest = cittest2;
        stat.A.second.tstudent(conta+3).tstat = statsttest2.tstat;
        stat.A.second.tstudent(conta+3).df = statsttest2.df;
        stat.A.second.tstudent(conta+3).sd = statsttest2.sd;
        stat.A.second.tstudent(conta+3).tstudent = veredict1;
        stat.A.second.tstudent(conta+3).sucess = success;
        stat.A.second.tstudent(conta+3).fail = fail;
        stat.A.second.tstudent(conta+3).data = (input_2yr(:,{'PatientCode','TWL','BMI'}));
        
        stat.A.second.Ftest(conta+3).fieldName = 'BMI';
        stat.A.second.Ftest(conta+3).fstat = stats_var.fstat;
        stat.A.second.Ftest(conta+3).df1 = stats_var.df1;
        stat.A.second.Ftest(conta+3).df2 = stats_var.df2;
        stat.A.second.Ftest(conta+3).h_ftest = h_var;
        stat.A.second.Ftest(conta+3).pval_ftest = p_var;
        stat.A.second.Ftest(conta+3).ci_ftest = ci_var;
        stat.A.second.Ftest(conta+3).ftest = veredict2;
        stat.A.second.Ftest(conta+3).sucess = success;
        stat.A.second.Ftest(conta+3).fail = fail;
        stat.A.second.Ftest(conta+3).data = (input_2yr(:,{'PatientCode','TWL','Age_'}));
        
        clearvars success fail

        % Dyslipidemia
        yes = input_2yr.dyslip == 1;
        no = input_2yr.dyslip == 0;

        figure('Name','Dyslipidemia 2yr')
        histogram(table2array(input_2yr(yes,{'TWL'})),...
            'BinMethod', 'integer', 'FaceColor', '#E19B00',...
            'EdgeColor','none')
        hold on
        histogram(table2array(input_2yr(no,{'TWL'})),...
            'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
        legend('yes','no')
        xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
            hold off
        if keep_fig == true
            title(['Distribution of success class and dyslipedemia (2yrs, n=',num2str(dsize),',',id_twl,')'])
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of success class and dyslipedemia (2yrs, n=',num2str(dsize),',',id_twl,')',fttype]);
        end
        
        tbl_cross = array2table(crosstab(yes,(input_2yr{:,'TWL'}))); % 1 = dyslipidemia
        %The rows in table correspond to the 2 distinct values in TWL, and
        %the columns correspond to the comorbidity.
        tbl_cross.Properties.RowNames = {'no_dyslip','dyslip'}; 
        tbl_cross.Properties.VariableNames = {'fail','success'};

        [h,p,stats] = fishertest(tbl_cross);
        if (h == 1)
            % Rejects null hypothesis
            veredict3 = ['Association dyslipidemia&TWL, ',...
                tbl_cross.Properties.RowNames{2}, ' have about ',...
                num2str(stats.OddsRatio),' greater odds of being ', ...
                tbl_cross.Properties.VariableNames{2}];
        else
            veredict3 = 'No association dyslipidemia&TWL';
        end

        stat.A.second.Fisher(conta+2).fieldName = 'Dyslipidemia by class';
        stat.A.second.Fisher(conta+2).h = h;
        stat.A.second.Fisher(conta+2).pval = p;
        stat.A.second.Fisher(conta+2).OddsRatio = stats.OddsRatio;
        stat.A.second.Fisher(conta+2).ci = stats.ConfidenceInterval;
        stat.A.second.Fisher(conta+2).Fisher = veredict3;
        stat.A.second.Fisher(conta+2).yes = yes;
        stat.A.second.Fisher(conta+2).no = no;
        stat.A.second.Fisher(conta+2).contingencyTable = tbl_cross;
        stat.A.second.Fisher(conta+2).data = (input_2yr(:,{'PatientCode','TWL','dyslip'}));
    
        clearvars yes no

        % Diabetes
        yes = input_2yr.diabet == 1;
        no = input_2yr.diabet == 0;

        figure('Name','Diabetes 2 yrs')
        histogram(table2array(input_2yr(yes,{'TWL'})),...
            'BinMethod', 'integer', 'FaceColor', '#E19B00',...
            'EdgeColor','none')
        hold on
        histogram(table2array(input_2yr(no,{'TWL'})),...
            'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
        legend('yes','no')
        xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
            hold off
        if keep_fig == true
            title(['Distribution of success class and diabetes (2yr, n=',...
                num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of success class and diabetes (2yr, n=',num2str(dsize),',',id_twl,')',fttype]);
        end
        
        tbl_cross = array2table(crosstab(yes,(input_2yr{:,'TWL'}))); % 1 = diabetes
        %The rows in table correspond to the 2 distinct values in TWL, and
        %the columns correspond to the comorbidity.
        tbl_cross.Properties.RowNames = {'no_diabetes','diabetes'}; 
        tbl_cross.Properties.VariableNames = {'fail','success'};

        [h,p,stats] = fishertest(tbl_cross);
        if (h == 1)
            % Rejects null hypothesis
            veredict3 = ['Association diabetes&TWL, ',...
                tbl_cross.Properties.RowNames{2}, ' have about ',...
                num2str(stats.OddsRatio),' greater odds of being ', ...
                tbl_cross.Properties.VariableNames{2}];
        else
            veredict3 = 'No association diabetes&TWL';
        end

        stat.A.second.Fisher(conta+3).fieldName = 'Diabetes by class';
        stat.A.second.Fisher(conta+3).h = h;
        stat.A.second.Fisher(conta+3).pval = p;
        stat.A.second.Fisher(conta+3).OddsRatio = stats.OddsRatio;
        stat.A.second.Fisher(conta+3).ci = stats.ConfidenceInterval;
        stat.A.second.Fisher(conta+3).Fisher = veredict3;
        stat.A.second.Fisher(conta+3).yes = yes;
        stat.A.second.Fisher(conta+3).no = no;
        stat.A.second.Fisher(conta+3).contingencyTable = tbl_cross;
        stat.A.second.Fisher(conta+3).data = (input_2yr(:,{'PatientCode','TWL','diabet'}));

        clearvars yes no

        % Hypertension
        yes = input_2yr.hypert == 1;
        no = input_2yr.hypert == 0;

        figure('Name','Hypertension 2yrs')
        histogram(table2array(input_2yr(yes,{'TWL'})),...
            'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
        histogram(table2array(input_2yr(no,{'TWL'})),...
            'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
        legend('yes','no')
        xlabel('Success class (1=succeed,0=failed)')
            hold off
        if keep_fig == true
            title(['Distribution of success class and hypertension (2yr, n=',...
                num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of success class and hypertension (2yr, n=',num2str(dsize),',',id_twl,')',fttype]);
        end   
        
        tbl_cross = array2table(crosstab(yes,(input_2yr{:,'TWL'}))); % 1 = hypertension
        %The rows in table correspond to the 2 distinct values in TWL, and
        %the columns correspond to the comorbidity.
        tbl_cross.Properties.RowNames = {'no_hypert','hypert'}; 
        tbl_cross.Properties.VariableNames = {'fail','success'};

        [h,p,stats] = fishertest(tbl_cross);
        if (h == 1)
            % Rejects null hypothesis
            veredict3 = ['Association hypertension&TWL, ',...
                tbl_cross.Properties.RowNames{2}, ' have about ',...
                num2str(stats.OddsRatio),' greater odds of being ', ...
                tbl_cross.Properties.VariableNames{2}];
        else
            veredict3 = 'No association hypertension&TWL';
        end

        stat.A.second.Fisher(conta+4).fieldName = 'Hypertension by class';
        stat.A.second.Fisher(conta+4).h = h;
        stat.A.second.Fisher(conta+4).pval = p;
        stat.A.second.Fisher(conta+4).OddsRatio = stats.OddsRatio;
        stat.A.second.Fisher(conta+4).ci = stats.ConfidenceInterval;
        stat.A.second.Fisher(conta+4).Fisher = veredict3;
        stat.A.second.Fisher(conta+4).yes = yes;
        stat.A.second.Fisher(conta+4).no = no;
        stat.A.second.Fisher(conta+4).contingencyTable = tbl_cross;
        stat.A.second.Fisher(conta+4).data = (input_2yr(:,{'PatientCode','TWL','hypert'}));

        clearvars yes no
    end
end

close all
%% Standarize Screening variables 1st yr

screen_std = {'rokenpy'};

for ppp=1:size(screen_std,2)
    var_screen = screen_std{ppp};
    s = (max(table2array(input_1yr(:,var_screen))) - min(table2array(input_1yr(:,var_screen))))/4; % sigma = (RefUp - RefLow)/4
    m = min(table2array(input_1yr(:,var_screen))) + 2*s;        % mu = RefLow + 2s
    for pp=1:size(input_1yr,1)
        input_1yr(pp,var_screen) = array2table((table2array(input_1yr(pp,var_screen))-m)/s);
    end
end
clearvars s m pp ppp var_screen

mat_1yr = table2array(input_1yr);

%% Identify variables with few measurements after 1 year

year = 1;
coordinates = x_SetnesKaymak2011(input_1yr,year,comp_foldername);
thr = height(input_1yr) - round(coordinates(1,1));
clear year

% Plots
f1=sum(~isnan(mat_1yr));

var_gradient = {'#ff0039','#fe0038','#fe0137','#fe0236','#fe0335',...
    '#fe0334','#fe0433','#fe0532','#fe0632','#fe0731','#fe0730', ...
    '#fe082f','#fe092e','#fe0a2d','#fe0a2c','#fe0b2c','#fe0c2b',...
    '#fe0d2a','#fe0e29','#fe0e28','#fe0f27','#fe1026','#fe1126',...
    '#fe1225','#fe1224','#fe1323','#fe1422','#fe1521','#fe1520',...
    '#fe1620','#fe171f','#fe181e','#fe191d','#fe191c','#fe1a1b',...
    '#fe1b1a','#fe1c19','#fe1d19','#fe1d18','#fe1e17','#fe1f16',...
    '#fe2015','#fe2014','#fe2113','#fe2213','#fe2312','#fe2411',...
    '#fe2410','#fe250f','#fe260e','#fe270d','#fe280d','#fe280c',...
    '#fe290b','#fe2a0a','#fe2b09','#fe2b08','#fe2c07','#fe2d07',...
    '#fe2e06','#fe2f05','#fd3302','#fe3301'};

var_gradient_rgb = zeros(length(var_gradient),3);
    
for i = 1:length(var_gradient)
    color_str = var_gradient{i};
    color_bar = sscanf(color_str(2:end),'%2x%2x%2x',[1 3])/255;
    var_gradient_rgb(i,:) = color_bar;
end

figure('Name','Features per patient gathered (1 year)');
subplot(2,2,[1,3])
h1 = bar(f1,'EdgeColor','none','FaceColor','flat');
h1.CData = var_gradient_rgb;
hline = refline([0 thr]);
hline.Color = 'w'; hline.LineWidth = 3;
xlabel(['Features, n=',num2str(length(f1))],'FontSize',16);
ylabel('Frequency/patient','FontSize',16);
title('Frequency before removing (1 yr)','FontSize',16);

idx_delete = find(sum(~isnan(mat_1yr)) < thr);
mat_1yr(:,idx_delete) = [];
input_1yr(:,idx_delete) = [];

% Result of previous exclusions
f2=sum(~isnan(mat_1yr));

var_gradient = {'#173f5f','#174060','#184261','#184462','#194663',...
    '#1a4865','#1a4a66','#1b4c67','#1c4e68','#1c5069','#1d526b',...
    '#1e546c','#1e566d','#1f586e','#205a70','#205c71','#215e72',...
    '#226073','#226274','#236476','#246677','#246878','#256a79',...
    '#266c7a','#266e7c','#27707d','#28727e','#28747f','#297681',...
    '#2a7882','#2a7a83','#2b7c84','#2c7e85','#2c8087','#2d8288',...
    '#2e8489','#2e868a','#2f888b','#308a8d','#308c8e','#318e8f',...
    '#329090','#329292','#339493','#349694','#349895','#359a96',...
    '#369c98','#369e99','#37a09a','#38a29b','#38a49c','#39a69e',...
    '#3aa89f','#3baba1','#3caea3'};

clearvars var_gradient_rgb

for i = 1:length(var_gradient)
    color_str = var_gradient{i};
    color_bar = sscanf(color_str(2:end),'%2x%2x%2x',[1 3])/255;
    var_gradient_rgb(i,:) = color_bar;
end

subplot(2,2,2)
h2 = bar(f2,'EdgeColor','none','FaceColor','flat');
h2.CData = var_gradient_rgb;
hline = refline([0 thr]);
hline.Color = 'w';hline.LineWidth = 3;
xlabel(['Features, n=',num2str(length(f2))],'FontSize',16);
ylabel('Frequency/patient','FontSize',16);
title('Frequency after removing (1 yr)','FontSize',16);

subplot(2,2,4)
h3=imagesc(isnan(mat_1yr(:,2:end)));
colormap('bone')
% colorbar
ylabel('White = missing, Black = no missing','FontSize',16)
xlabel('Variables','FontSize',16)
title('Remaining missing entries','FontSize',16);
if keep_fig == true
    fig = gcf; ax = gca;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 28; %Width
    fig.Position(4) = 30; %Height
    ax.FontSize = 14;
   saveas(gcf,[comp_foldername '/Features 1 yr, n=',num2str(length(f2)),fttype]);
end 
clearvars f2 f1 hline d

close all

%% Standarize Screening variables 2nd yr

for ppp=1:size(screen_std,2)
    var_screen = screen_std{ppp};
    s = (max(table2array(input_2yr(:,var_screen))) - min(table2array(input_2yr(:,var_screen))))/4; % sigma = (RefUp - RefLow)/4
    m = min(table2array(input_2yr(:,var_screen))) + 2*s;        % mu = RefLow + 2s
    for pp=1:size(input_2yr,1)
        input_2yr(pp,var_screen) = array2table((table2array(input_2yr(pp,var_screen))-m)/s);
    end
end
clearvars s m pp ppp var_screen

mat_2yr = table2array(input_2yr);

%% Identify variables with few measurements after 2 years

year = 2;
coordinates2 = x_SetnesKaymak2011(input_2yr,year,comp_foldername);
thr2 = height(input_2yr) - round(coordinates2(1,1));
clearvars year

var_gradient = {'#ff0039','#fe0038','#fe0137','#fe0236','#fe0335',...
    '#fe0434','#fe0433','#fe0532','#fe0631','#fe0731','#fe0830',...
    '#fe082f','#fe092e','#fe0a2d','#fe0b2c','#fe0c2b','#fe0c2a',...
    '#fe0d29','#fe0e29','#fe0f28','#fe1027','#fe1126','#fe1125',...
    '#fe1224','#fe1323','#fe1422','#fe1521','#fe1521','#fe1620',...
    '#fe171f','#fe181e','#fe191d','#fe191c','#fe1a1b','#fe1b1a',...
    '#fe1c19','#fe1d19','#fe1d18','#fe1e17','#fe1f16','#fe2015',...
    '#fe2114','#fe2213','#fe2212','#fe2311','#fe2411','#fe2510',...
    '#fe260f','#fe260e','#fe270d','#fe280c','#fe290b','#fe2a0a',...
    '#fe2a09','#fe2b09','#fe2c08','#fe2d07','#fe2e06','#fe2e05',...
    '#fe2f04','#fe3003','#fe3102','#fe3301'};

var_gradient_rgb = zeros(length(var_gradient),3);
    
for i = 1:length(var_gradient)
    color_str = var_gradient{i};
    color_bar = sscanf(color_str(2:end),'%2x%2x%2x',[1 3])/255;
    var_gradient_rgb(i,:) = color_bar;
end

% Plot
f3=sum(~isnan(mat_2yr));
figure('Name','Features per patient gathered (2 years)');
subplot(2,2,[1,3])
h4 = bar(f3,'EdgeColor','none','FaceColor','flat'); %#ok<*NASGU>
h4.CData = var_gradient_rgb;
hline = refline([0 thr2]);
hline.Color = 'w'; hline.LineWidth = 3;
xlabel(['Features, n=',num2str(length(f3))],'FontSize',16);
ylabel('Frequency/patient','FontSize',16);
title('Frequency before removing (2 yrs)','FontSize',16);

idx_delete = find(sum(~isnan(mat_2yr)) < thr2);
mat_2yr(:,idx_delete) = [];
input_2yr(:,idx_delete) = [];

var_gradient = {'#173f5f','#174060','#184261','#184462','#194663',...
    '#1a4865','#1a4a66','#1b4c67','#1c4e68','#1c5069','#1d526b',...
    '#1e546c','#1e566d','#1f586e','#205a70','#205c71','#215e72',...
    '#226073','#226274','#236476','#246677','#246878','#256a79',...
    '#266c7a','#266e7c','#27707d','#28727e','#28747f','#297681',...
    '#2a7882','#2a7a83','#2b7c84','#2c7e85','#2c8087','#2d8288',...
    '#2e8489','#2e868a','#2f888b','#308a8d','#308c8e','#318e8f',...
    '#329090','#329292','#339493','#349694','#349895','#359a96',...
    '#369c98','#369e99','#37a09a','#38a29b','#38a49c','#39a69e',...
    '#3aa89f','#3aaaa0','#3caea3'};

var_gradient_rgb = zeros(length(var_gradient),3);
    
for i = 1:length(var_gradient)
    color_str = var_gradient{i};
    color_bar = sscanf(color_str(2:end),'%2x%2x%2x',[1 3])/255;
    var_gradient_rgb(i,:) = color_bar;
end


% Result of previous exclusions
f4=sum(~isnan(mat_2yr));

subplot(2,2,2)
h5 = bar(f4,'EdgeColor','none','FaceColor','flat');
h5.CData = var_gradient_rgb;
hline = refline([0 thr2]);
hline.Color = 'w'; hline.LineWidth = 3;
xlabel(['Features, n=',num2str(length(f4))],'FontSize',16);
ylabel('Frequency/patient','FontSize',16);
title('Frequency after removing (2 yrs)','FontSize',16);

subplot(2,2,4)
h6=imagesc(isnan(mat_2yr(:,2:end)));
colormap('bone')
% colorbar
ylabel('White = missing, Black = no missing','FontSize',16)
xlabel('Variables','FontSize',16)
title('Remaining missing entries','FontSize',16);
if keep_fig == true
    fig = gcf; ax = gca;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 28; %Width
    fig.Position(4) = 30; %Height
    ax.FontSize = 14;
   saveas(gcf,[comp_foldername '/Features 2 yr, n=',num2str(length(f4)),fttype]);
end

clearvars f3 f4 hline

close all

%% Delete-wise method
% Criterion D (1st part) - All patients

% 1 year
tbl_in_1yr = input_1yr(~any(ismissing(input_1yr,NaN),2),:);

if hist_dist == true
    dsize = size(tbl_in_1yr,1);
    % Gender
    male = tbl_in_1yr.Sex == 1;
    female = tbl_in_1yr.Sex == 0;
    
    figure('Name','Gender 1yr - critD')
    histogram(tbl_in_1yr{:,'TWL'}(female),'BinMethod', 'integer',...
        'FaceColor', '#CF0047','EdgeColor','none','FaceAlpha',0.4);
        hold on
    histogram(tbl_in_1yr{:,'TWL'}(male),'BinMethod', 'integer',...
        'FaceColor', '#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    legend('female','male','FontSize',12)
    xlabel('Output variable (%TWL)','FontSize',16,...
        'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of TWL among gender (1yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among gender (1yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    figure('Name','Age by sex - crit D')
        hold on
    histogram(tbl_in_1yr{:,'Age_'}(female),'BinMethod', 'integer',...
        'FaceColor','#CF0047','FaceAlpha',0.4,'EdgeColor','none');
    histogram(tbl_in_1yr{:,'Age_'}(male),'BinMethod', 'integer',...
        'FaceColor','#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    hold off
    legend('female','male','FontSize',12)
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
    if keep_fig == true
        title(['Distribution of ages by sex (1yr, n=',...
        num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of ages (1 yr n=', num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    figure('Name','Normplot Age - crit C')
    h1 = normplot(tbl_in_1yr{:,'Age_'});
    h1(1).LineWidth = 2;h1(2).LineWidth = 2;h1(3).LineWidth = 2;
    h1(1).MarkerEdgeColor = '#E5CC00';
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age (1 year, n=', ...
        num2str(dsize),',',id_twl,'C)',fttype]);
    
    figure('Name','Normplots Age by sex')
    hold on
    hm = normplot(tbl_in_1yr{:,'Age_'}(male));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_in_1yr{:,'Age_'}(female));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Male','','','Female'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age by sex (1 yr, n=', ...
        num2str(dsize),',',id_twl,'C)',fttype]);
    
    % Null hypothesis: The mean strengths for the two populations are equal.
    % A p-value less than the significance level indicates that you 
    % can reject the null hypothesis. In other words, the sample provides 
    % sufficient evidence to conclude that the population means are 
    % different. Below is the output for the analysis.
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_in_1yr{:,'Age_'}(male)),...
    (tbl_in_1yr{:,'Age_'}(female)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_in_1yr{:,'Age_'}(male)),...
    (tbl_in_1yr{:,'Age_'}(female)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.B.first.tstudent(conta+1).fieldName = 'Age by sex';
    stat.B.first.tstudent(conta+1).h_ttest = httest2;
    stat.B.first.tstudent(conta+1).pval_ttest = pttest2;
    stat.B.first.tstudent(conta+1).ci_ttest = cittest2;
    stat.B.first.tstudent(conta+1).tstat = statsttest2.tstat;
    stat.B.first.tstudent(conta+1).df = statsttest2.df;
    stat.B.first.tstudent(conta+1).sd = statsttest2.sd;
    stat.B.first.tstudent(conta+1).tstudent = veredict1;
    stat.B.first.tstudent(conta+1).male = male;
    stat.B.first.tstudent(conta+1).female = female;
    stat.B.first.tstudent(conta+1).data = (tbl_in_1yr(:,{'PatientCode','Sex','Age_'}));
    
    
    stat.B.first.Ftest(conta+1).fieldName = 'Age by sex';
    stat.B.first.Ftest(conta+1).fstat = stats_var.fstat;
    stat.B.first.Ftest(conta+1).df1 = stats_var.df1;
    stat.B.first.Ftest(conta+1).df2 = stats_var.df2;
    stat.B.first.Ftest(conta+1).h_ftest = h_var;
    stat.B.first.Ftest(conta+1).pval_ftest = p_var;
    stat.B.first.Ftest(conta+1).ci_ftest = ci_var;
    stat.B.first.Ftest(conta+1).ftest = veredict2;
    stat.B.first.Ftest(conta+1).male = male;
    stat.B.first.Ftest(conta+1).female = female;
    stat.B.first.Ftest(conta+1).data = (tbl_in_1yr(:,{'PatientCode','Sex','Age_'}));
    
    tbl_cross = array2table(crosstab(female,(tbl_in_1yr{:,'TWL'}))); % 1 = female
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the 2 distinct GENDERS.
    tbl_cross.Properties.RowNames = {'male','female'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association gender & age ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association gender & age';
    end
    
    stat.B.first.Fisher(conta+1).fieldName = 'Sex by class';
    stat.B.first.Fisher(conta+1).h = h;
    stat.B.first.Fisher(conta+1).pval = p;
    stat.B.first.Fisher(conta+1).OddsRatio = stats.OddsRatio;
    stat.B.first.Fisher(conta+1).ci = stats.ConfidenceInterval;
    stat.B.first.Fisher(conta+1).Fisher = veredict3;
    stat.B.first.Fisher(conta+1).contingencyTable = tbl_cross;
    stat.B.first.Fisher(conta+1).yes = 'TWL';
    stat.B.first.Fisher(conta+1).no = 'TWL';
    stat.B.first.Fisher(conta+1).data = (tbl_in_1yr(:,{'PatientCode','TWL','Age_'}));
    
    clearvars male female
    
     % BMI
    success = tbl_in_1yr.TWL == 1;
    fail = tbl_in_1yr.TWL == 0;
    
    figure('Name','BMI 1st yr - critC')
    histogram(table2array(tbl_in_1yr(success,{'BMI'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(tbl_in_1yr(fail,{'BMI'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('BMI (kg/m^2)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class among BMI (1yr, n=', num2str(dsize),',',id_twl,'C)'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among BMI (1yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    figure('Name','Normplots BMI')
    hold on
    hm = normplot(tbl_in_1yr{:,'BMI'}(success));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_in_1yr{:,'BMI'}(fail));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Success','','','Fail'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of BMI by class (1 yr, n=', ...
        num2str(dsize),',',id_twl,'C)',fttype]);
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_in_1yr{:,'BMI'}(success)),...
    (tbl_in_1yr{:,'BMI'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_in_1yr{:,'BMI'}(success)),...
    (tbl_in_1yr{:,'BMI'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.B.first.tstudent(conta+2).fieldName = 'BMI';
    stat.B.first.tstudent(conta+2).h_ttest = httest2;
    stat.B.first.tstudent(conta+2).pval_ttest = pttest2;
    stat.B.first.tstudent(conta+2).ci_ttest = cittest2;
    stat.B.first.tstudent(conta+2).tstat = statsttest2.tstat;
    stat.B.first.tstudent(conta+2).df = statsttest2.df;
    stat.B.first.tstudent(conta+2).sd = statsttest2.sd;
    stat.B.first.tstudent(conta+2).tstudent = veredict1;
    stat.B.first.tstudent(conta+2).sucess = success;
    stat.B.first.tstudent(conta+2).fail = fail;
    stat.B.first.tstudent(conta+2).data = (tbl_in_1yr(:,{'PatientCode','TWL','BMI'}));
    
    [pWhitney,hWhitney,statsWhitney] = ranksum((tbl_in_1yr{:,'BMI'}(success)),...
    (tbl_in_1yr{:,'BMI'}(fail)));
    if (hWhitney == 1)
        % rejects the null hypothesis
        veredictW = 'Diff sig (median)'; 
    else
        veredictW = 'No diff sig (median)';
    end
    
    stat.B.first.Whitney(conta+1).fieldName = 'BMI';
    stat.B.first.Whitney(conta+1).h = hWhitney;
    stat.B.first.Whitney(conta+1).p = pWhitney;
    stat.B.first.Whitney(conta+1).ranksum = statsWhitney.ranksum;
    stat.B.first.Whitney(conta+1).Whitney = veredictW;
    stat.B.first.Whitney(conta+1).sucess = success;
    stat.B.first.Whitney(conta+1).fail = fail;
    stat.B.first.Whitney(conta+1).data = (input_1yr(:,{'PatientCode','TWL','BMI'}));
    
    stat.B.first.Ftest(conta+2).fieldName = 'BMI';
    stat.B.first.Ftest(conta+2).fstat = stats_var.fstat;
    stat.B.first.Ftest(conta+2).df1 = stats_var.df1;
    stat.B.first.Ftest(conta+2).df2 = stats_var.df2;
    stat.B.first.Ftest(conta+2).h_ftest = h_var;
    stat.B.first.Ftest(conta+2).pval_ftest = p_var;
    stat.B.first.Ftest(conta+2).ci_ftest = ci_var;
    stat.B.first.Ftest(conta+2).ftest = veredict2;
    stat.B.first.Ftest(conta+2).sucess = success;
    stat.B.first.Ftest(conta+2).fail = fail;
    stat.B.first.Ftest(conta+2).data = (tbl_in_1yr(:,{'PatientCode','TWL','BMI'}));
    
    % Age
    figure('Name','Age 1 yr - critC')
    histogram(table2array(tbl_in_1yr(success,{'Age_'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(tbl_in_1yr(fail,{'Age_'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
        hold off
        if keep_fig == true
            title(['Distribution of ages among classes (1yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of ages among classes (1yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
        end
    
    figure('Name','Normplot Age (class)')
    hold on
    hm = normplot(tbl_in_1yr{:,'Age_'}(success));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_in_1yr{:,'Age_'}(fail));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Success','','','Fail'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of TWL by sex (1 yr, n=', ...
        num2str(dsize),',',id_twl,'C)',fttype]);
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_in_1yr{:,'Age_'}(fail)),...
    (tbl_in_1yr{:,'Age_'}(success)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_in_1yr{:,'Age_'}(fail)),...
    (tbl_in_1yr{:,'Age_'}(success)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.B.first.tstudent(conta+3).fieldName = 'Age by class';
    stat.B.first.tstudent(conta+3).h_ttest = httest2;
    stat.B.first.tstudent(conta+3).pval_ttest = pttest2;
    stat.B.first.tstudent(conta+3).ci_ttest = cittest2;
    stat.B.first.tstudent(conta+3).tstat = statsttest2.tstat;
    stat.B.first.tstudent(conta+3).df = statsttest2.df;
    stat.B.first.tstudent(conta+3).sd = statsttest2.sd;
    stat.B.first.tstudent(conta+3).tstudent = veredict1;
    stat.B.first.tstudent(conta+3).sucess = success;
    stat.B.first.tstudent(conta+3).fail = fail;
    stat.B.first.tstudent(conta+3).data = (tbl_in_1yr(:,{'PatientCode','Sex','Age_'}));
    
    stat.B.first.Ftest(conta+3).fieldName = 'Age by class';
    stat.B.first.Ftest(conta+3).fstat = stats_var.fstat;
    stat.B.first.Ftest(conta+3).df1 = stats_var.df1;
    stat.B.first.Ftest(conta+3).df2 = stats_var.df2;
    stat.B.first.Ftest(conta+3).h_ftest = h_var;
    stat.B.first.Ftest(conta+3).pval_ftest = p_var;
    stat.B.first.Ftest(conta+3).ci_ftest = ci_var;
    stat.B.first.Ftest(conta+3).ftest = veredict2;
    stat.B.first.Ftest(conta+3).success = success;
    stat.B.first.Ftest(conta+3).fail = fail;
    stat.B.first.Ftest(conta+3).data = (tbl_in_1yr(:,{'PatientCode','Sex','Age_'}));
    
    clearvars success fail
    
    % Dyslipidemia
    yes = tbl_in_1yr.dyslip == max(tbl_in_1yr.dyslip);
    no = tbl_in_1yr.dyslip == min(tbl_in_1yr.dyslip);
    
    figure('Name','Dislipidemia 1 yr')
    histogram(table2array(tbl_in_1yr(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_in_1yr(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and dyslipedemia (1yr, n=', num2str(dsize),',',id_twl,'C)'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and dyslipedemia (1yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_in_1yr{:,'TWL'}))); % 1 = dyslipidemia
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_dyslip','dyslip'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association dyslipidemia&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association dyslipidemia&TWL';
    end
   
    stat.B.first.Fisher(conta+2).fieldName = 'Dyslipidemia by class';
    stat.B.first.Fisher(conta+2).h = h;
    stat.B.first.Fisher(conta+2).pval = p;
    stat.B.first.Fisher(conta+2).OddsRatio = stats.OddsRatio;
    stat.B.first.Fisher(conta+2).ci = stats.ConfidenceInterval;
    stat.B.first.Fisher(conta+2).Fisher = veredict3;
    stat.B.first.Fisher(conta+2).yes = yes;
    stat.B.first.Fisher(conta+2).no = no;
    stat.B.first.Fisher(conta+2).contingencyTable = tbl_cross;
    stat.B.first.Fisher(conta+2).data = (tbl_in_1yr(:,{'PatientCode','TWL','dyslip'}));
    
    clearvars yes no
    
    % Diabetes
    yes = tbl_in_1yr.diabet == max(tbl_in_1yr.diabet);
    no = tbl_in_1yr.diabet == min(tbl_in_1yr.diabet);
    
    figure('Name','Diabetes 1 yr')
    histogram(table2array(tbl_in_1yr(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_in_1yr(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
     if keep_fig == true
        title(['Distribution of success class and diabetes (1yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and diabetes (1yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
     end
     
    tbl_cross = array2table(crosstab(yes,(tbl_in_1yr{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_diabetes','diabetes'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association diabetes&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association diabetes&TWL';
    end
    
    stat.B.first.Fisher(conta+3).fieldName = 'Diabetes by class';
    stat.B.first.Fisher(conta+3).h = h;
    stat.B.first.Fisher(conta+3).pval = p;
    stat.B.first.Fisher(conta+3).OddsRatio = stats.OddsRatio;
    stat.B.first.Fisher(conta+3).ci = stats.ConfidenceInterval;
    stat.B.first.Fisher(conta+3).Fisher = veredict3;
    stat.B.first.Fisher(conta+3).yes = yes;
    stat.B.first.Fisher(conta+3).no = no;
    stat.B.first.Fisher(conta+3).contingencyTable = tbl_cross;
    stat.B.first.Fisher(conta+3).data = (tbl_in_1yr(:,{'PatientCode','TWL','diabet'}));
    
    clearvars yes no
    
    % Hypertension
    yes = tbl_in_1yr.hypert == max(tbl_in_1yr.hypert);
    no = tbl_in_1yr.hypert == min(tbl_in_1yr.hypert);
    
    figure('Name','Hypertension 1 yrs')
    histogram(table2array(tbl_in_1yr(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_in_1yr(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and hypertension (1yr, n=',num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and hypertension (1yr, n=',num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_in_1yr{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_hypert','hypert'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association hypertension&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association hypertension&TWL';
    end
    
    stat.B.first.Fisher(conta+4).fieldName = 'Hypertension by class';
    stat.B.first.Fisher(conta+4).h = h;
    stat.B.first.Fisher(conta+4).pval = p;
    stat.B.first.Fisher(conta+4).OddsRatio = stats.OddsRatio;
    stat.B.first.Fisher(conta+4).ci = stats.ConfidenceInterval;
    stat.B.first.Fisher(conta+4).Fisher = veredict3;
    stat.B.first.Fisher(conta+4).yes = yes;
    stat.B.first.Fisher(conta+4).no = no;
    stat.B.first.Fisher(conta+4).contingencyTable = tbl_cross;
    stat.B.first.Fisher(conta+4).data = (tbl_in_1yr(:,{'PatientCode','TWL','hypert'}));
    
    clearvars yes no
end

% 2 years
tbl_in_2yr = input_2yr(~any(ismissing(input_2yr,NaN),2),:);

if hist_dist == true
    dsize = size(tbl_in_2yr,1);
    % Gender
    male = tbl_in_2yr.Sex == 1;
    female = tbl_in_2yr.Sex == 0;
    
    figure('Name','Gender 2yr - critC')
    histogram(tbl_in_2yr{:,'TWL'}(female),'BinMethod', 'integer',...
        'FaceColor', '#CF0047','EdgeColor','none','FaceAlpha',0.4);
        hold on
    histogram(tbl_in_2yr{:,'TWL'}(male),'BinMethod', 'integer',...
        'FaceColor', '#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    legend('female','male','FontSize',12)
    xlabel('Output variable (%TWL)','FontSize',16,...
        'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of TWL among gender (2yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among gender (2yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    figure('Name','Age by sex 2yr - critC')
        hold on
    histogram(tbl_in_2yr{:,'Age_'}(female),'BinMethod', 'integer',...
        'FaceColor','#CF0047','FaceAlpha',0.4,'EdgeColor','none');
    histogram(tbl_in_2yr{:,'Age_'}(male),'BinMethod', 'integer',...
        'FaceColor','#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    hold off
    legend('female','male','FontSize',12)
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
    if keep_fig == true
        title(['Distribution of ages by sex (2yr, n=',...
        num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of ages (2 yr n=', num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    figure('Name','Normplot Age')
    h1 = normplot(tbl_in_2yr{:,'Age_'});
    h1(1).LineWidth = 2;h1(2).LineWidth = 2;h1(3).LineWidth = 2;
    h1(1).MarkerEdgeColor = '#E5CC00';
    xlabel('Age at surgery (2 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age (2 year, n=', ...
        num2str(dsize),',',id_twl,'C)',fttype]);

    figure('Name','Normplots Age by sex')
    hold on
    hm = normplot(tbl_in_2yr{:,'Age_'}(male));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_in_2yr{:,'Age_'}(female));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Male','','','Female'},'Location','southeast')
    xlabel('Age at surgery (2 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age by sex (2 yr, n=', ...
        num2str(dsize),',',id_twl,'C)',fttype]);
    
    % Null hypothesis: The mean strengths for the two populations are equal.
    % A p-value less than the significance level indicates that you 
    % can reject the null hypothesis. In other words, the sample provides 
    % sufficient evidence to conclude that the population means are 
    % different. Below is the output for the analysis.
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_in_2yr{:,'Age_'}(male)),...
    (tbl_in_2yr{:,'Age_'}(female)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_in_2yr{:,'Age_'}(male)),...
    (tbl_in_2yr{:,'Age_'}(female)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.B.second.tstudent(conta+1).fieldName = 'Age by sex';
    stat.B.second.tstudent(conta+1).h_ttest = httest2;
    stat.B.second.tstudent(conta+1).pval_ttest = pttest2;
    stat.B.second.tstudent(conta+1).ci_ttest = cittest2;
    stat.B.second.tstudent(conta+1).tstat = statsttest2.tstat;
    stat.B.second.tstudent(conta+1).df = statsttest2.df;
    stat.B.second.tstudent(conta+1).sd = statsttest2.sd;
    stat.B.second.tstudent(conta+1).tstudent = veredict1;
    stat.B.second.tstudent(conta+1).male = male;
    stat.B.second.tstudent(conta+1).female = female;
    stat.B.second.tstudent(conta+1).data = (tbl_in_2yr(:,{'PatientCode','Sex','Age_'}));

    stat.B.second.Ftest(conta+1).fieldName = 'Age by sex';
    stat.B.second.Ftest(conta+1).fstat = stats_var.fstat;
    stat.B.second.Ftest(conta+1).df1 = stats_var.df1;
    stat.B.second.Ftest(conta+1).df2 = stats_var.df2;
    stat.B.second.Ftest(conta+1).h_ftest = h_var;
    stat.B.second.Ftest(conta+1).pval_ftest = p_var;
    stat.B.second.Ftest(conta+1).ci_ftest = ci_var;
    stat.B.second.Ftest(conta+1).ftest = veredict2;
    stat.B.second.Ftest(conta+1).male = male;
    stat.B.second.Ftest(conta+1).female = female;
    stat.B.second.Ftest(conta+1).data = (tbl_in_2yr(:,{'PatientCode','Sex','Age_'}));

    tbl_cross = array2table(crosstab(female,(tbl_in_2yr{:,'TWL'}))); % 1 = female
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the 2 distinct GENDERS.
    tbl_cross.Properties.RowNames = {'male','female'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};

    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association gender & age ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association gender & age';
    end

    stat.B.second.Fisher(conta+1).fieldName = 'Sex by class';
    stat.B.second.Fisher(conta+1).h = h;
    stat.B.second.Fisher(conta+1).pval = p;
    stat.B.second.Fisher(conta+1).OddsRatio = stats.OddsRatio;
    stat.B.second.Fisher(conta+1).ci = stats.ConfidenceInterval;
    stat.B.second.Fisher(conta+1).Fisher = veredict3;
    stat.B.second.Fisher(conta+1).contingencyTable = tbl_cross;
    stat.B.second.Fisher(conta+1).yes = 'TWL';
    stat.B.second.Fisher(conta+1).no = 'TWL';
    stat.B.second.Fisher(conta+1).data = (tbl_in_2yr(:,{'PatientCode','TWL','Age_'}));

    clearvars male female
    
     % BMI
    success = tbl_in_2yr.TWL == 1;
    fail = tbl_in_2yr.TWL == 0;
    
    figure('Name','BMI 1st yr - critC')
    histogram(table2array(tbl_in_2yr(success,{'BMI'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(tbl_in_2yr(fail,{'BMI'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('BMI (kg/m^2)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class among BMI min (2yr, n=', num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among BMI min (2yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    figure('Name','Normplots BMI')
    hold on
    hm = normplot(tbl_in_2yr{:,'BMI'}(success));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_in_2yr{:,'BMI'}(fail));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Success','','','Fail'},'Location','southeast')
    xlabel('Age at surgery (2 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of BMI by class (2 yr, n=', ...
        num2str(dsize),',',id_twl,'C)',fttype]);

    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_in_2yr{:,'BMI'}(success)),...
    (tbl_in_2yr{:,'BMI'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end

    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_in_2yr{:,'BMI'}(success)),...
    (tbl_in_2yr{:,'BMI'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.B.second.tstudent(conta+2).fieldName = 'BMI';
    stat.B.second.tstudent(conta+2).h_ttest = httest2;
    stat.B.second.tstudent(conta+2).pval_ttest = pttest2;
    stat.B.second.tstudent(conta+2).ci_ttest = cittest2;
    stat.B.second.tstudent(conta+2).tstat = statsttest2.tstat;
    stat.B.second.tstudent(conta+2).df = statsttest2.df;
    stat.B.second.tstudent(conta+2).sd = statsttest2.sd;
    stat.B.second.tstudent(conta+2).tstudent = veredict1;
    stat.B.second.tstudent(conta+2).sucess = success;
    stat.B.second.tstudent(conta+2).fail = fail;
    stat.B.second.tstudent(conta+2).data = (tbl_in_2yr(:,{'PatientCode','TWL','BMI'}));
    
    [pWhitney,hWhitney,statsWhitney] = ranksum((tbl_in_2yr{:,'BMI'}(success)),...
    (tbl_in_2yr{:,'BMI'}(fail)));
    if (hWhitney == 1)
        % rejects the null hypothesis
        veredictW = 'Diff sig (median)'; 
    else
        veredictW = 'No diff sig (median)';
    end
    
    stat.B.first.Whitney(conta+1).fieldName = 'BMI';
    stat.B.first.Whitney(conta+1).h = hWhitney;
    stat.B.first.Whitney(conta+1).p = pWhitney;
    stat.B.first.Whitney(conta+1).ranksum = statsWhitney.ranksum;
    stat.B.first.Whitney(conta+1).Whitney = veredictW;
    stat.B.first.Whitney(conta+1).sucess = success;
    stat.B.first.Whitney(conta+1).fail = fail;
    stat.B.first.Whitney(conta+1).data = (tbl_in_2yr(:,{'PatientCode','TWL','BMI'}));
    
    stat.B.second.Ftest(conta+2).fieldName = 'BMI';
    stat.B.second.Ftest(conta+2).fstat = stats_var.fstat;
    stat.B.second.Ftest(conta+2).df1 = stats_var.df1;
    stat.B.second.Ftest(conta+2).df2 = stats_var.df2;
    stat.B.second.Ftest(conta+2).h_ftest = h_var;
    stat.B.second.Ftest(conta+2).pval_ftest = p_var;
    stat.B.second.Ftest(conta+2).ci_ftest = ci_var;
    stat.B.second.Ftest(conta+2).ftest = veredict2;
    stat.B.second.Ftest(conta+2).sucess = success;
    stat.B.second.Ftest(conta+2).fail = fail;
    stat.B.second.Ftest(conta+2).data = (tbl_in_2yr(:,{'PatientCode','TWL','BMI'}));

    % Age
    figure('Name','Age 1 yr - critC')
    histogram(table2array(tbl_in_2yr(success,{'Age_'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(tbl_in_2yr(fail,{'Age_'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
        hold off
        if keep_fig == true
            title(['Distribution of ages among classes (2yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of ages among classes (2yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
        end
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_in_2yr{:,'Age_'}(fail)),...
    (tbl_in_2yr{:,'Age_'}(success)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_in_2yr{:,'Age_'}(fail)),...
    (tbl_in_2yr{:,'Age_'}(success)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.B.second.tstudent(conta+3).fieldName = 'Age by class';
    stat.B.second.tstudent(conta+3).h_ttest = httest2;
    stat.B.second.tstudent(conta+3).pval_ttest = pttest2;
    stat.B.second.tstudent(conta+3).ci_ttest = cittest2;
    stat.B.second.tstudent(conta+3).tstat = statsttest2.tstat;
    stat.B.second.tstudent(conta+3).df = statsttest2.df;
    stat.B.second.tstudent(conta+3).sd = statsttest2.sd;
    stat.B.second.tstudent(conta+3).tstudent = veredict1;
    stat.B.second.tstudent(conta+3).sucess = success;
    stat.B.second.tstudent(conta+3).fail = fail;
    stat.B.second.tstudent(conta+3).data = (tbl_in_2yr(:,{'PatientCode','Sex','Age_'}));
    
    stat.B.second.Ftest(conta+3).fieldName = 'Age by class';
    stat.B.second.Ftest(conta+3).fstat = stats_var.fstat;
    stat.B.second.Ftest(conta+3).df1 = stats_var.df1;
    stat.B.second.Ftest(conta+3).df2 = stats_var.df2;
    stat.B.second.Ftest(conta+3).h_ftest = h_var;
    stat.B.second.Ftest(conta+3).pval_ftest = p_var;
    stat.B.second.Ftest(conta+3).ci_ftest = ci_var;
    stat.B.second.Ftest(conta+3).ftest = veredict2;
    stat.B.second.Ftest(conta+3).success = success;
    stat.B.second.Ftest(conta+3).success = fail;
    stat.B.second.Ftest(conta+3).data = (tbl_in_2yr(:,{'PatientCode','Sex','Age_'}));
    
    clearvars success fail
    
    % Dyslipidemia
    yes = tbl_in_2yr.dyslip == max(tbl_in_2yr.dyslip);
    no = tbl_in_2yr.dyslip == min(tbl_in_2yr.dyslip);
    
    figure('Name','Dyslipidemia 1 yr')
    histogram(table2array(tbl_in_2yr(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_in_2yr(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and dyslipedemia (2yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and dyslipedemia (2yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_in_2yr{:,'TWL'}))); % 1 = dyslipidemia
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_dyslip','dyslip'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};

    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association dyslipidemia&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association dyslipidemia&TWL';
    end

    stat.B.second.Fisher(conta+2).fieldName = 'Dyslipidemia by class';
    stat.B.second.Fisher(conta+2).h = h;
    stat.B.second.Fisher(conta+2).pval = p;
    stat.B.second.Fisher(conta+2).OddsRatio = stats.OddsRatio;
    stat.B.second.Fisher(conta+2).ci = stats.ConfidenceInterval;
    stat.B.second.Fisher(conta+2).Fisher = veredict3;
    stat.B.second.Fisher(conta+2).yes = yes;
    stat.B.second.Fisher(conta+2).no = no;
    stat.B.second.Fisher(conta+2).contingencyTable = tbl_cross;
    stat.B.second.Fisher(conta+2).data = (tbl_in_2yr(:,{'PatientCode','TWL','dyslip'}));

    clearvars yes no
    
    % Diabetes
    yes = tbl_in_2yr.diabet == max(tbl_in_2yr.diabet);
    no = tbl_in_2yr.diabet == min(tbl_in_2yr.diabet);
    
    figure('Name','Diabetes 1 yr')
    histogram(table2array(tbl_in_2yr(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_in_2yr(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
    hold off
     if keep_fig == true
        title(['Distribution of success class and diabetes (2yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and diabetes (2yr, n=', num2str(dsize),',',id_twl,'C)',fttype]);
     end
    
    tbl_cross = array2table(crosstab(yes,(tbl_in_2yr{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_diabetes','diabetes'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};

    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association diabetes&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association diabetes&TWL';
    end

    stat.B.second.Fisher(conta+3).fieldName = 'Diabetes by class';
    stat.B.second.Fisher(conta+3).h = h;
    stat.B.second.Fisher(conta+3).pval = p;
    stat.B.second.Fisher(conta+3).OddsRatio = stats.OddsRatio;
    stat.B.second.Fisher(conta+3).ci = stats.ConfidenceInterval;
    stat.B.second.Fisher(conta+3).Fisher = veredict3;
    stat.B.second.Fisher(conta+3).yes = yes;
    stat.B.second.Fisher(conta+3).no = no;
    stat.B.second.Fisher(conta+3).contingencyTable = tbl_cross;
    stat.B.second.Fisher(conta+3).data = (tbl_in_2yr(:,{'PatientCode','TWL','diabet'}));

    clearvars yes no
    
    % Hypertension
    yes = tbl_in_2yr.hypert == max(tbl_in_2yr.hypert);
    no = tbl_in_2yr.hypert == min(tbl_in_2yr.hypert);
    
    figure('Name','Hypertension 1 yrs')
    histogram(table2array(tbl_in_2yr(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_in_2yr(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and hypertension (2yr, n=',num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and hypertension (2yr, n=',num2str(dsize),',',id_twl,'C)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_in_2yr{:,'TWL'}))); % 1 = hypertension
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_hypert','hypert'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};

    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association hypertension&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association hypertension&TWL';
    end

    stat.B.second.Fisher(conta+4).fieldName = 'Hypertension by class';
    stat.B.second.Fisher(conta+4).h = h;
    stat.B.second.Fisher(conta+4).pval = p;
    stat.B.second.Fisher(conta+4).OddsRatio = stats.OddsRatio;
    stat.B.second.Fisher(conta+4).ci = stats.ConfidenceInterval;
    stat.B.second.Fisher(conta+4).Fisher = veredict3;
    stat.B.second.Fisher(conta+4).yes = yes;
    stat.B.second.Fisher(conta+4).no = no;
    stat.B.second.Fisher(conta+4).contingencyTable = tbl_cross;
    stat.B.second.Fisher(conta+4).data = (tbl_in_2yr(:,{'PatientCode','TWL','hypert'}));

    clearvars yes no
end

close all

%% Get patients that went trough metabolic surgery and those who underwent bariatric
% Criterion D (2nd part)

% 1st year
idx_bariatric_1yr = (table2array(tbl_in_1yr(:,'Metabolic')))>=40.0;
tbl_1yr_bariatric = tbl_in_1yr(idx_bariatric_1yr,:);
tbl_1yr_metabolic = tbl_in_1yr(~idx_bariatric_1yr,:);
tbl_1yr_metabolic(:,{'Metabolic'}) = [];
tbl_1yr_bariatric(:,{'Metabolic'}) = [];
tbl_in_1yr(:,{'Metabolic'}) = [];

if hist_dist == true
    dsize = size(tbl_1yr_metabolic,1);
    % Gender
    male = tbl_1yr_metabolic.Sex == 1;
    female = tbl_1yr_metabolic.Sex == 0;
    
    figure('Name','Gender 1yr - critD')
    histogram(tbl_1yr_metabolic{:,'TWL'}(female),'BinMethod', 'integer',...
        'FaceColor', '#CF0047','EdgeColor','none','FaceAlpha',0.4);
    hold on
    histogram(tbl_1yr_metabolic{:,'TWL'}(male),'BinMethod', 'integer',...
        'FaceColor', '#0097CD','FaceAlpha',0.4,'EdgeColor','none');
%     [f1,xi1] = ksdensity(tbl_1yr_metabolic{:,'TWL'}(female));
%     [f2,xi2] = ksdensity(tbl_1yr_metabolic{:,'TWL'}(male));
%     p = plot(xi1,f1,xi2,f2);
%     p(1).Color = '#CF0047'; p(2).Color = '#0097CD';
%     p(1).LineWidth = 2; p(2).LineWidth = 2;
    legend('female','male','FontSize',12)
    xlabel('Output variable (%TWL)','FontSize',16,...
        'FontWeight','bold')
    hold off
    if keep_fig == true
        title(['Distribution of TWL among gender (1yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among gender (1yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    figure('Name','Age by sex 1yr - critD')
%     hold on
%     histogram(tbl_1yr_metabolic{:,'Age_'}(female),'BinMethod', 'integer',...
%         'FaceColor','#CF0047','FaceAlpha',0.4,'EdgeColor','none');
%     histogram(tbl_1yr_metabolic{:,'Age_'}(male),'BinMethod', 'integer',...
%         'FaceColor','#0097CD','FaceAlpha',0.4,'EdgeColor','none');
%     hold off
    [f1,xi1] = ksdensity(tbl_1yr_metabolic{:,'Age_'}(female));
    [f2,xi2] = ksdensity(tbl_1yr_metabolic{:,'Age_'}(male));
    p = plot(xi1,f1,xi2,f2);
    p(1).Color = '#CF0047'; p(2).Color = '#0097CD';
    p(1).LineWidth = 2; p(2).LineWidth = 2;
    legend('female','male','FontSize',12)
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
    if keep_fig == true
        title(['Distribution of ages by sex (1yr, n=',...
        num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of ages (1 yr n=', num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    figure('Name','Normplot Age - crit C')
    h1 = normplot(tbl_1yr_metabolic{:,'Age_'});
    h1(1).LineWidth = 2;h1(2).LineWidth = 2;h1(3).LineWidth = 2;
    h1(1).MarkerEdgeColor = '#E5CC00';
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age (1 year, n=', ...
        num2str(dsize),',',id_twl,'D)',fttype]);
    
    figure('Name','Normplots Age by sex')
    hold on
    hm = normplot(tbl_1yr_metabolic{:,'Age_'}(male));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_1yr_metabolic{:,'Age_'}(female));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Male','','','Female'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age by sex (1 yr, n=', ...
        num2str(dsize),',',id_twl,'D)',fttype]);
    
    % Null hypothesis: The mean strengths for the two populations are equal.
    % A p-value less than the significance level indicates that you 
    % can reject the null hypothesis. In other words, the sample provides 
    % sufficient evidence to conclude that the population means are 
    % different. Below is the output for the analysis.
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_1yr_metabolic{:,'Age_'}(male)),...
    (tbl_1yr_metabolic{:,'Age_'}(female)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_1yr_metabolic{:,'Age_'}(male)),...
    (tbl_1yr_metabolic{:,'Age_'}(female)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.C.first.tstudent(conta+1).fieldName = 'Age by sex';
    stat.C.first.tstudent(conta+1).h_ttest = httest2;
    stat.C.first.tstudent(conta+1).pval_ttest = pttest2;
    stat.C.first.tstudent(conta+1).ci_ttest = cittest2;
    stat.C.first.tstudent(conta+1).tstat = statsttest2.tstat;
    stat.C.first.tstudent(conta+1).df = statsttest2.df;
    stat.C.first.tstudent(conta+1).sd = statsttest2.sd;
    stat.C.first.tstudent(conta+1).tstudent = veredict1;
    stat.C.first.tstudent(conta+1).male = male;
    stat.C.first.tstudent(conta+1).female = female;
    stat.C.first.tstudent(conta+1).data = (tbl_1yr_metabolic(:,{'PatientCode','Sex','Age_'}));
    
    stat.C.first.Ftest(conta+1).fieldName = 'Age by sex';
    stat.C.first.Ftest(conta+1).fstat = stats_var.fstat;
    stat.C.first.Ftest(conta+1).df1 = stats_var.df1;
    stat.C.first.Ftest(conta+1).df2 = stats_var.df2;
    stat.C.first.Ftest(conta+1).h_ftest = h_var;
    stat.C.first.Ftest(conta+1).pval_ftest = p_var;
    stat.C.first.Ftest(conta+1).ci_ftest = ci_var;
    stat.C.first.Ftest(conta+1).ftest = veredict2;
    stat.C.first.Ftest(conta+1).male = male;
    stat.C.first.Ftest(conta+1).female = female;
    stat.C.first.Ftest(conta+1).data = (tbl_1yr_metabolic(:,{'PatientCode','Sex','Age_'}));
    
    tbl_cross = array2table(crosstab(female,(tbl_1yr_metabolic{:,'TWL'}))); % 1 = female
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the 2 distinct GENDERS.
    tbl_cross.Properties.RowNames = {'male','female'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association gender & age ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association gender & age';
    end
    
    stat.C.first.Fisher(conta+1).fieldName = 'Gender by class';
    stat.C.first.Fisher(conta+1).h = h;
    stat.C.first.Fisher(conta+1).pval = p;
    stat.C.first.Fisher(conta+1).OddsRatio = stats.OddsRatio;
    stat.C.first.Fisher(conta+1).ci = stats.ConfidenceInterval;
    stat.C.first.Fisher(conta+1).Fisher = veredict3;
    stat.C.first.Fisher(conta+1).contingencyTable = tbl_cross;
    stat.C.first.Fisher(conta+1).yes = 'TWL';
    stat.C.first.Fisher(conta+1).no = 'TWL';
    stat.C.first.Fisher(conta+1).data = (tbl_1yr_metabolic(:,{'PatientCode','TWL','Age_'}));
    
    clearvars male female
    
     % BMI
    success = tbl_1yr_metabolic.TWL == 1;
    fail = tbl_1yr_metabolic.TWL == 0;
    
    figure('Name','BMI 1st yr - critD')
%     histogram(table2array(tbl_1yr_metabolic(success,{'BMI'})),...
%         'BinMethod', 'integer', 'FaceColor', '#E19B00',...
%         'EdgeColor','none')
%         hold on
%     histogram(table2array(tbl_1yr_metabolic(fail,{'BMI'})),...
%         'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
%         'EdgeColor','none');
    [f1,xi1] = ksdensity(table2array(tbl_1yr_metabolic(success,{'BMI'})));
    [f2,xi2] = ksdensity(table2array(tbl_1yr_metabolic(fail,{'BMI'})));
    p = plot(xi1,f1,xi2,f2);
    p(1).Color = '#E19B00'; p(2).Color = '#CF0047';
    p(1).LineWidth = 2; p(2).LineWidth = 2;
    legend('success','failure')
    xlabel('BMI (kg/m^2)','FontSize',16,'FontWeight','bold')
%         hold off
    if keep_fig == true
        title(['Distribution of success class among BMI min (1yr, n=', num2str(dsize),',',id_twl,')'],'FontSize',18)
        ylabel('Density','FontSize',16,'FontWeight','bold')
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among BMI min (1yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    figure('Name','Normplots BMI')
    hold on
    hm = normplot(tbl_1yr_metabolic{:,'BMI'}(success));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_1yr_metabolic{:,'BMI'}(fail));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Success','','','Fail'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of BMI by class (1 yr, n=', ...
        num2str(dsize),',',id_twl,'D)',fttype]);
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_1yr_metabolic{:,'BMI'}(success)),...
    (tbl_1yr_metabolic{:,'BMI'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_1yr_metabolic{:,'BMI'}(success)),...
    (tbl_1yr_metabolic{:,'BMI'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.C.first.tstudent(conta+2).fieldName = 'BMI';
    stat.C.first.tstudent(conta+2).h_ttest = httest2;
    stat.C.first.tstudent(conta+2).pval_ttest = pttest2;
    stat.C.first.tstudent(conta+2).ci_ttest = cittest2;
    stat.C.first.tstudent(conta+2).tstat = statsttest2.tstat;
    stat.C.first.tstudent(conta+2).df = statsttest2.df;
    stat.C.first.tstudent(conta+2).sd = statsttest2.sd;
    stat.C.first.tstudent(conta+2).tstudent = veredict1;
    stat.C.first.tstudent(conta+2).sucess = success;
    stat.C.first.tstudent(conta+2).fail = fail;
    stat.C.first.tstudent(conta+2).data = (tbl_1yr_metabolic(:,{'PatientCode','TWL','BMI'}));
    
    [pWhitney,hWhitney,statsWhitney] = ranksum((tbl_1yr_metabolic{:,'BMI'}(success)),...
    (tbl_1yr_metabolic{:,'BMI'}(fail)));
    if (hWhitney == 1)
        % rejects the null hypothesis
        veredictW = 'Diff sig (median)'; 
    else
        veredictW = 'No diff sig (median)';
    end
    
    stat.C.first.Whitney(conta+1).fieldName = 'BMI';
    stat.C.first.Whitney(conta+1).h = hWhitney;
    stat.C.first.Whitney(conta+1).p = pWhitney;
    stat.C.first.Whitney(conta+1).ranksum = statsWhitney.ranksum;
    stat.C.first.Whitney(conta+1).Whitney = veredictW;
    stat.C.first.Whitney(conta+1).sucess = success;
    stat.C.first.Whitney(conta+1).fail = fail;
    stat.C.first.Whitney(conta+1).data = (tbl_1yr_metabolic(:,{'PatientCode','TWL','BMI'}));
    
    stat.C.first.Ftest(conta+2).fieldName = 'BMI';
    stat.C.first.Ftest(conta+2).fstat = stats_var.fstat;
    stat.C.first.Ftest(conta+2).df1 = stats_var.df1;
    stat.C.first.Ftest(conta+2).df2 = stats_var.df2;
    stat.C.first.Ftest(conta+2).h_ftest = h_var;
    stat.C.first.Ftest(conta+2).pval_ftest = p_var;
    stat.C.first.Ftest(conta+2).ci_ftest = ci_var;
    stat.C.first.Ftest(conta+2).ftest = veredict2;
    stat.C.first.Ftest(conta+2).sucess = success;
    stat.C.first.Ftest(conta+2).fail = fail;
    stat.C.first.Ftest(conta+2).data = (tbl_1yr_metabolic(:,{'PatientCode','TWL','BMI'}));
    
    % Age
    figure('Name','Age 1 yr - critC')
%     histogram(table2array(tbl_1yr_metabolic(success,{'Age_'})),...
%         'BinMethod', 'integer', 'FaceColor', '#E19B00',...
%         'EdgeColor','none')
%         hold on
%     histogram(table2array(tbl_1yr_metabolic(fail,{'Age_'})),...
%         'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
%         'EdgeColor','none');
    [f1,xi1] = ksdensity(table2array(tbl_1yr_metabolic(success,{'Age_'})));
    [f2,xi2] = ksdensity(table2array(tbl_1yr_metabolic(fail,{'Age_'})));
    p = plot(xi1,f1,xi2,f2);
    p(1).Color = '#E19B00'; p(2).Color = '#CF0047';
    p(1).LineWidth = 2; p(2).LineWidth = 2;
    legend('success','failure')
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
%         hold off
    if keep_fig == true
        title(['Distribution of ages among classes (1yr, n=',...
        num2str(dsize),',',id_twl,')'],'FontSize',18)
        ylabel('Density','FontSize',16,'FontWeight','bold')
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of ages among classes (1yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_1yr_metabolic{:,'Age_'}(success)),...
    (tbl_1yr_metabolic{:,'Age_'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_1yr_metabolic{:,'Age_'}(success)),...
    (tbl_1yr_metabolic{:,'Age_'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.C.first.tstudent(conta+3).fieldName = 'Age and Class';
    stat.C.first.tstudent(conta+3).h_ttest = httest2;
    stat.C.first.tstudent(conta+3).pval_ttest = pttest2;
    stat.C.first.tstudent(conta+3).ci_ttest = cittest2;
    stat.C.first.tstudent(conta+3).tstat = statsttest2.tstat;
    stat.C.first.tstudent(conta+3).df = statsttest2.df;
    stat.C.first.tstudent(conta+3).sd = statsttest2.sd;
    stat.C.first.tstudent(conta+3).tstudent = veredict1;
    stat.C.first.tstudent(conta+3).sucess = success;
    stat.C.first.tstudent(conta+3).fail = fail;
    stat.C.first.tstudent(conta+3).data = (input_1yr(:,{'PatientCode','TWL','BMI'}));
    
    stat.C.first.Ftest(conta+3).fieldName = 'BMI';
    stat.C.first.Ftest(conta+3).fstat = stats_var.fstat;
    stat.C.first.Ftest(conta+3).df1 = stats_var.df1;
    stat.C.first.Ftest(conta+3).df2 = stats_var.df2;
    stat.C.first.Ftest(conta+3).h_ftest = h_var;
    stat.C.first.Ftest(conta+3).pval_ftest = p_var;
    stat.C.first.Ftest(conta+3).ci_ftest = ci_var;
    stat.C.first.Ftest(conta+3).ftest = veredict2;
    stat.C.first.Ftest(conta+3).sucess = success;
    stat.C.first.Ftest(conta+3).fail = fail;
    stat.C.first.Ftest(conta+3).data = (input_1yr(:,{'PatientCode','TWL','Age_'}));
    
    clearvars success fail
    
    % Dyslipidemia
    yes = tbl_1yr_metabolic.dyslip == max(tbl_1yr_metabolic.dyslip);
    no = tbl_1yr_metabolic.dyslip == min(tbl_1yr_metabolic.dyslip);
    
    figure('Name','Dyslipidemia 1 yr')
    histogram(table2array(tbl_1yr_metabolic(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_1yr_metabolic(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and dyslipedemia (1yr, n=', num2str(dsize),',',id_twl,'D)'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and dyslipedemia (1yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_1yr_metabolic{:,'TWL'}))); % 1 = dyslipidemia
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_dyslip','dyslip'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association dyslipidemia&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association dyslipidemia&TWL';
    end
   
    stat.C.first.Fisher(conta+2).fieldName = 'Dyslipidemia by class';
    stat.C.first.Fisher(conta+2).h = h;
    stat.C.first.Fisher(conta+2).pval = p;
    stat.C.first.Fisher(conta+2).OddsRatio = stats.OddsRatio;
    stat.C.first.Fisher(conta+2).ci = stats.ConfidenceInterval;
    stat.C.first.Fisher(conta+2).Fisher = veredict3;
    stat.C.first.Fisher(conta+2).yes = yes;
    stat.C.first.Fisher(conta+2).no = no;
    stat.C.first.Fisher(conta+2).contingencyTable = tbl_cross;
    stat.C.first.Fisher(conta+2).data = (tbl_1yr_metabolic(:,{'PatientCode','TWL','dyslip'}));
    
    clearvars yes no
    
    % Diabetes
    yes = tbl_1yr_metabolic.diabet == max(tbl_1yr_metabolic.diabet);
    no = tbl_1yr_metabolic.diabet == min(tbl_1yr_metabolic.diabet);
    
    figure('Name','Diabetes 1 yr')
    histogram(table2array(tbl_1yr_metabolic(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_1yr_metabolic(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
     if keep_fig == true
        title(['Distribution of success class and diabetes (1yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and diabetes (1yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
     end
     
    tbl_cross = array2table(crosstab(yes,(tbl_1yr_metabolic{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_diabetes','diabetes'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association diabetes&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association diabetes&TWL';
    end
    
    stat.C.first.Fisher(conta+3).fieldName = 'Diabetes by class';
    stat.C.first.Fisher(conta+3).h = h;
    stat.C.first.Fisher(conta+3).pval = p;
    stat.C.first.Fisher(conta+3).OddsRatio = stats.OddsRatio;
    stat.C.first.Fisher(conta+3).ci = stats.ConfidenceInterval;
    stat.C.first.Fisher(conta+3).Fisher = veredict3;
    stat.C.first.Fisher(conta+3).yes = yes;
    stat.C.first.Fisher(conta+3).no = no;
    stat.C.first.Fisher(conta+3).contingencyTable = tbl_cross;
    stat.C.first.Fisher(conta+3).data = (tbl_1yr_metabolic(:,{'PatientCode','TWL','diabet'}));
    
    clearvars yes no
    
    % Hypertension
    yes = tbl_1yr_metabolic.hypert == max(tbl_1yr_metabolic.hypert);
    no = tbl_1yr_metabolic.hypert == min(tbl_1yr_metabolic.hypert);
    
    figure('Name','Hypertension 1 yrs')
    histogram(table2array(tbl_1yr_metabolic(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_1yr_metabolic(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and hypertension (1yr, n=',num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and hypertension (1yr, n=',num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_1yr_metabolic{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_hypert','hypert'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association hypertension&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association hypertension&TWL';
    end
    
    stat.C.first.Fisher(conta+4).fieldName = 'Hypertension by class';
    stat.C.first.Fisher(conta+4).h = h;
    stat.C.first.Fisher(conta+4).pval = p;
    stat.C.first.Fisher(conta+4).OddsRatio = stats.OddsRatio;
    stat.C.first.Fisher(conta+4).ci = stats.ConfidenceInterval;
    stat.C.first.Fisher(conta+4).Fisher = veredict3;
    stat.C.first.Fisher(conta+4).yes = yes;
    stat.C.first.Fisher(conta+4).no = no;
    stat.C.first.Fisher(conta+4).contingencyTable = tbl_cross;
    stat.C.first.Fisher(conta+4).data = (tbl_1yr_metabolic(:,{'PatientCode','TWL','hypert'}));
    
    clearvars yes no
end

close all

if hist_dist == true
    dsize = size(tbl_1yr_bariatric,1);
    % Gender
    male = tbl_1yr_bariatric.Sex == 1;
    female = tbl_1yr_bariatric.Sex == 0;
    
    figure('Name','Gender 1yr - critD')
    histogram(tbl_1yr_bariatric{:,'TWL'}(female),'BinMethod', 'integer',...
        'FaceColor', '#CF0047','EdgeColor','none','FaceAlpha',0.4);
        hold on
    histogram(tbl_1yr_bariatric{:,'TWL'}(male),'BinMethod', 'integer',...
        'FaceColor', '#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    legend('female','male','FontSize',12)
    xlabel('Output variable (%TWL)','FontSize',16,...
        'FontWeight','bold')
        hold off
    if keep_fig == true
        ylabel('Density','FontSize',16,'FontWeight','bold')
        title(['Distribution of TWL among gender (1yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among gender (1yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    figure('Name','Age by sex 1yr - critD')
%     hold on
%     histogram(tbl_1yr_bariatric{:,'Age_'}(female),'BinMethod', 'integer',...
%         'FaceColor','#CF0047','FaceAlpha',0.4,'EdgeColor','none');
%     histogram(tbl_1yr_bariatric{:,'Age_'}(male),'BinMethod', 'integer',...
%         'FaceColor','#0097CD','FaceAlpha',0.4,'EdgeColor','none');
%     hold off
    [f1,xi1] = ksdensity(tbl_1yr_bariatric{:,'Age_'}(female));
    [f2,xi2] = ksdensity(tbl_1yr_bariatric{:,'Age_'}(male));
    p = plot(xi1,f1,xi2,f2);
    p(1).Color = '#CF0047'; p(2).Color = '#0097CD';
    p(1).LineWidth = 2; p(2).LineWidth = 2;
    legend('female','male','FontSize',12)
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
    if keep_fig == true
        title(['Distribution of ages by sex (1yr, n=',...
        num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of ages (1 yr n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    figure('Name','Normplot Age - crit C')
    h1 = normplot(tbl_1yr_bariatric{:,'Age_'});
    h1(1).LineWidth = 2;h1(2).LineWidth = 2;h1(3).LineWidth = 2;
    h1(1).MarkerEdgeColor = '#E5CC00';
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age (1 year, n=', ...
        num2str(dsize),',',id_twl,'Dbar)',fttype]);
    
    figure('Name','Normplots Age by sex')
    hold on
    hm = normplot(tbl_1yr_bariatric{:,'Age_'}(male));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_1yr_bariatric{:,'Age_'}(female));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Male','','','Female'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age by sex (1 yr, n=', ...
        num2str(dsize),',',id_twl,'Dbar)',fttype]);
    
    % Null hypothesis: The mean strengths for the two populations are equal.
    % A p-value less than the significance level indicates that you 
    % can reject the null hypothesis. In other words, the sample provides 
    % sufficient evidence to conclude that the population means are 
    % different. Below is the output for the analysis.
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_1yr_bariatric{:,'Age_'}(male)),...
    (tbl_1yr_bariatric{:,'Age_'}(female)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_1yr_bariatric{:,'Age_'}(male)),...
    (tbl_1yr_bariatric{:,'Age_'}(female)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.D.first.tstudent(conta+1).fieldName = 'Age by sex';
    stat.D.first.tstudent(conta+1).h_ttest = httest2;
    stat.D.first.tstudent(conta+1).pval_ttest = pttest2;
    stat.D.first.tstudent(conta+1).ci_ttest = cittest2;
    stat.D.first.tstudent(conta+1).tstat = statsttest2.tstat;
    stat.D.first.tstudent(conta+1).df = statsttest2.df;
    stat.D.first.tstudent(conta+1).sd = statsttest2.sd;
    stat.D.first.tstudent(conta+1).tstudent = veredict1;
    stat.D.first.tstudent(conta+1).male = male;
    stat.D.first.tstudent(conta+1).female = female;
    stat.D.first.tstudent(conta+1).data = (tbl_1yr_bariatric(:,{'PatientCode','Sex','Age_'}));
    
    stat.D.first.Ftest(conta+1).fieldName = 'Age by sex';
    stat.D.first.Ftest(conta+1).fstat = stats_var.fstat;
    stat.D.first.Ftest(conta+1).df1 = stats_var.df1;
    stat.D.first.Ftest(conta+1).df2 = stats_var.df2;
    stat.D.first.Ftest(conta+1).h_ftest = h_var;
    stat.D.first.Ftest(conta+1).pval_ftest = p_var;
    stat.D.first.Ftest(conta+1).ci_ftest = ci_var;
    stat.D.first.Ftest(conta+1).ftest = veredict2;
    stat.D.first.Ftest(conta+1).male = male;
    stat.D.first.Ftest(conta+1).female = female;
    stat.D.first.Ftest(conta+1).data = (tbl_1yr_bariatric(:,{'PatientCode','Sex','Age_'}));
    
    tbl_cross = array2table(crosstab(female,(tbl_1yr_bariatric{:,'TWL'}))); % 1 = female
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the 2 distinct GENDERS.
    tbl_cross.Properties.RowNames = {'male','female'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association gender & TWL ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association gender & TWL';
    end
    
    stat.D.first.Fisher(conta+1).fieldName = 'Gender by class';
    stat.D.first.Fisher(conta+1).h = h;
    stat.D.first.Fisher(conta+1).pval = p;
    stat.D.first.Fisher(conta+1).OddsRatio = stats.OddsRatio;
    stat.D.first.Fisher(conta+1).ci = stats.ConfidenceInterval;
    stat.D.first.Fisher(conta+1).Fisher = veredict3;
    stat.D.first.Fisher(conta+1).contingencyTable = tbl_cross;
    stat.D.first.Fisher(conta+1).yes = 'TWL';
    stat.D.first.Fisher(conta+1).no = 'TWL';
    stat.D.first.Fisher(conta+1).data = (tbl_1yr_bariatric(:,{'PatientCode','TWL','Age_'}));
    
    clearvars male female
    
     % BMI
    success = tbl_1yr_bariatric.TWL == 1;
    fail = tbl_1yr_bariatric.TWL == 0;
    
    figure('Name','BMI 1st yr - critD')
    histogram(table2array(tbl_1yr_bariatric(success,{'BMI'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(tbl_1yr_bariatric(fail,{'BMI'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('BMI (kg/m^2)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class among BMI min (1yr, n=', num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among BMI min (1yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    figure('Name','Normplots BMI')
    hold on
    hm = normplot(tbl_1yr_bariatric{:,'BMI'}(success));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_1yr_bariatric{:,'BMI'}(fail));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Success','','','Fail'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of BMI by class (1 yr, n=', ...
        num2str(dsize),',',id_twl,'Dbar)',fttype]);
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_1yr_bariatric{:,'BMI'}(success)),...
    (tbl_1yr_bariatric{:,'BMI'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_1yr_bariatric{:,'BMI'}(success)),...
    (tbl_1yr_bariatric{:,'BMI'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.D.first.tstudent(conta+2).fieldName = 'BMI';
    stat.D.first.tstudent(conta+2).h_ttest = httest2;
    stat.D.first.tstudent(conta+2).pval_ttest = pttest2;
    stat.D.first.tstudent(conta+2).ci_ttest = cittest2;
    stat.D.first.tstudent(conta+2).tstat = statsttest2.tstat;
    stat.D.first.tstudent(conta+2).df = statsttest2.df;
    stat.D.first.tstudent(conta+2).sd = statsttest2.sd;
    stat.D.first.tstudent(conta+2).tstudent = veredict1;
    stat.D.first.tstudent(conta+2).sucess = success;
    stat.D.first.tstudent(conta+2).fail = fail;
    stat.D.first.tstudent(conta+2).data = (tbl_1yr_bariatric(:,{'PatientCode','TWL','BMI'}));
    
    [pWhitney,hWhitney,statsWhitney] = ranksum((tbl_1yr_bariatric{:,'BMI'}(success)),...
    (tbl_1yr_bariatric{:,'BMI'}(fail)));
    if (hWhitney == 1)
        % rejects the null hypothesis
        veredictW = 'Diff sig (median)'; 
    else
        veredictW = 'No diff sig (median)';
    end
    
    stat.D.first.Whitney(conta+1).fieldName = 'BMI';
    stat.D.first.Whitney(conta+1).h = hWhitney;
    stat.D.first.Whitney(conta+1).p = pWhitney;
    stat.D.first.Whitney(conta+1).ranksum = statsWhitney.ranksum;
    stat.D.first.Whitney(conta+1).Whitney = veredictW;
    stat.D.first.Whitney(conta+1).sucess = success;
    stat.D.first.Whitney(conta+1).fail = fail;
    stat.D.first.Whitney(conta+1).data = (tbl_1yr_bariatric(:,{'PatientCode','TWL','BMI'}));
    
    stat.D.first.Ftest(conta+2).fieldName = 'BMI';
    stat.D.first.Ftest(conta+2).fstat = stats_var.fstat;
    stat.D.first.Ftest(conta+2).df1 = stats_var.df1;
    stat.D.first.Ftest(conta+2).df2 = stats_var.df2;
    stat.D.first.Ftest(conta+2).h_ftest = h_var;
    stat.D.first.Ftest(conta+2).pval_ftest = p_var;
    stat.D.first.Ftest(conta+2).ci_ftest = ci_var;
    stat.D.first.Ftest(conta+2).ftest = veredict2;
    stat.D.first.Ftest(conta+2).sucess = success;
    stat.D.first.Ftest(conta+2).fail = fail;
    stat.D.first.Ftest(conta+2).data = (tbl_1yr_bariatric(:,{'PatientCode','TWL','BMI'}));
    
    % Age
    figure('Name','Age 1 yr - critD')
    histogram(table2array(tbl_1yr_bariatric(success,{'Age_'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(tbl_1yr_bariatric(fail,{'Age_'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
        hold off
        if keep_fig == true
            title(['Distribution of ages among classes (1yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of ages among classes (1yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
        end
    
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_1yr_bariatric{:,'Age_'}(success)),...
    (tbl_1yr_bariatric{:,'Age_'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_1yr_bariatric{:,'Age_'}(success)),...
    (tbl_1yr_bariatric{:,'Age_'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.D.first.tstudent(conta+3).fieldName = 'Age and class';
    stat.D.first.tstudent(conta+3).h_ttest = httest2;
    stat.D.first.tstudent(conta+3).pval_ttest = pttest2;
    stat.D.first.tstudent(conta+3).ci_ttest = cittest2;
    stat.D.first.tstudent(conta+3).tstat = statsttest2.tstat;
    stat.D.first.tstudent(conta+3).df = statsttest2.df;
    stat.D.first.tstudent(conta+3).sd = statsttest2.sd;
    stat.D.first.tstudent(conta+3).tstudent = veredict1;
    stat.D.first.tstudent(conta+3).sucess = success;
    stat.D.first.tstudent(conta+3).fail = fail;
    stat.D.first.tstudent(conta+3).data = (tbl_1yr_bariatric(:,{'PatientCode','TWL','BMI'}));
    
    stat.D.first.Ftest(conta+3).fieldName = 'Age and class';
    stat.D.first.Ftest(conta+3).fstat = stats_var.fstat;
    stat.D.first.Ftest(conta+3).df1 = stats_var.df1;
    stat.D.first.Ftest(conta+3).df2 = stats_var.df2;
    stat.D.first.Ftest(conta+3).h_ftest = h_var;
    stat.D.first.Ftest(conta+3).pval_ftest = p_var;
    stat.D.first.Ftest(conta+3).ci_ftest = ci_var;
    stat.D.first.Ftest(conta+3).ftest = veredict2;
    stat.D.first.Ftest(conta+3).sucess = success;
    stat.D.first.Ftest(conta+3).fail = fail;
    stat.D.first.Ftest(conta+3).data = (tbl_1yr_bariatric(:,{'PatientCode','TWL','Age_'}));
        
    clearvars success fail
    
    % Dyslipidemia
    yes = tbl_1yr_bariatric.dyslip == max(tbl_1yr_bariatric.dyslip);
    no = tbl_1yr_bariatric.dyslip == min(tbl_1yr_bariatric.dyslip);
    
    figure('Name','Dislipidemia 1 yr')
    histogram(table2array(tbl_1yr_bariatric(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_1yr_bariatric(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and dyslipedemia (1yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and dyslipedemia (1yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_1yr_bariatric{:,'TWL'}))); % 1 = dyslipidemia
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_dyslip','dyslip'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association dyslipidemia&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association dyslipidemia&TWL';
    end
   
    stat.D.first.Fisher(conta+2).fieldName = 'Dyslipidemia by class';
    stat.D.first.Fisher(conta+2).h = h;
    stat.D.first.Fisher(conta+2).pval = p;
    stat.D.first.Fisher(conta+2).OddsRatio = stats.OddsRatio;
    stat.D.first.Fisher(conta+2).ci = stats.ConfidenceInterval;
    stat.D.first.Fisher(conta+2).Fisher = veredict3;
    stat.D.first.Fisher(conta+2).yes = yes;
    stat.D.first.Fisher(conta+2).no = no;
    stat.D.first.Fisher(conta+2).contingencyTable = tbl_cross;
    stat.D.first.Fisher(conta+2).data = (tbl_1yr_bariatric(:,{'PatientCode','TWL','dyslip'}));
    
    clearvars yes no
    
    % Diabetes
    yes = tbl_1yr_bariatric.diabet == max(tbl_1yr_bariatric.diabet);
    no = tbl_1yr_bariatric.diabet == min(tbl_1yr_bariatric.diabet);
    
    figure('Name','Diabetes 1 yr')
    histogram(table2array(tbl_1yr_bariatric(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_1yr_bariatric(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
     if keep_fig == true
        title(['Distribution of success class and diabetes (1yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and diabetes (1yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
     end
    
    tbl_cross = array2table(crosstab(yes,(tbl_1yr_bariatric{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_diabetes','diabetes'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association diabetes&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association diabetes&TWL';
    end
    
    stat.D.first.Fisher(conta+3).fieldName = 'Diabetes by class';
    stat.D.first.Fisher(conta+3).h = h;
    stat.D.first.Fisher(conta+3).pval = p;
    stat.D.first.Fisher(conta+3).OddsRatio = stats.OddsRatio;
    stat.D.first.Fisher(conta+3).ci = stats.ConfidenceInterval;
    stat.D.first.Fisher(conta+3).Fisher = veredict3;
    stat.D.first.Fisher(conta+3).yes = yes;
    stat.D.first.Fisher(conta+3).no = no;
    stat.D.first.Fisher(conta+3).contingencyTable = tbl_cross;
    stat.D.first.Fisher(conta+3).data = (tbl_1yr_bariatric(:,{'PatientCode','TWL','diabet'}));
    
    clearvars yes no
    
    % Hypertension
    yes = tbl_1yr_bariatric.hypert == max(tbl_1yr_bariatric.hypert);
    no = tbl_1yr_bariatric.hypert == min(tbl_1yr_bariatric.hypert);
    
    figure('Name','Hypertension 1 yrs')
    histogram(table2array(tbl_1yr_bariatric(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_1yr_bariatric(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and hypertension (1yr, n=',num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and hypertension (1yr, n=',num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_1yr_bariatric{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_hypert','hypert'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association hypertension&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association hypertension&TWL';
    end
    
    stat.D.first.Fisher(conta+4).fieldName = 'Hypertension by class';
    stat.D.first.Fisher(conta+4).h = h;
    stat.D.first.Fisher(conta+4).pval = p;
    stat.D.first.Fisher(conta+4).OddsRatio = stats.OddsRatio;
    stat.D.first.Fisher(conta+4).ci = stats.ConfidenceInterval;
    stat.D.first.Fisher(conta+4).Fisher = veredict3;
    stat.D.first.Fisher(conta+4).yes = yes;
    stat.D.first.Fisher(conta+4).no = no;
    stat.D.first.Fisher(conta+4).contingencyTable = tbl_cross;
    stat.D.first.Fisher(conta+4).data = (tbl_1yr_bariatric(:,{'PatientCode','TWL','hypert'}));
    
    clearvars yes no
end

close all

% 2nd year

idx_bariatric_2yr = (table2array(tbl_in_2yr(:,'Metabolic')))>=40.0;
tbl_2yr_bariatric = tbl_in_2yr(idx_bariatric_2yr,:);
tbl_2yr_metabolic = tbl_in_2yr(~idx_bariatric_2yr,:);
tbl_2yr_bariatric(:,{'Metabolic'}) = [];
tbl_2yr_metabolic(:,{'Metabolic'}) = [];
tbl_in_2yr(:,{'Metabolic'}) = [];

if hist_dist == true
    dsize = size(tbl_2yr_metabolic,1);
    % Gender
    male = tbl_2yr_metabolic.Sex == 1;
    female = tbl_2yr_metabolic.Sex == 0;
    
    figure('Name','Gender 1yr - critD')
    histogram(tbl_2yr_metabolic{:,'TWL'}(female),'BinMethod', 'integer',...
        'FaceColor', '#CF0047','EdgeColor','none','FaceAlpha',0.4);
        hold on
    histogram(tbl_2yr_metabolic{:,'TWL'}(male),'BinMethod', 'integer',...
        'FaceColor', '#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    legend('female','male','FontSize',12)
    xlabel('Output variable (%TWL)','FontSize',16,...
        'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of TWL among gender (2yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among gender (2yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    figure('Name','Age by sex 2yr - critD')
%     hold on
%     histogram(tbl_2yr_metabolic{:,'Age_'}(female),'BinMethod', 'integer',...
%         'FaceColor','#CF0047','FaceAlpha',0.4,'EdgeColor','none');
%     histogram(tbl_2yr_metabolic{:,'Age_'}(male),'BinMethod', 'integer',...
%         'FaceColor','#0097CD','FaceAlpha',0.4,'EdgeColor','none');
%     hold off
    [f1,xi1] = ksdensity(tbl_2yr_metabolic{:,'Age_'}(female));
    [f2,xi2] = ksdensity(tbl_2yr_metabolic{:,'Age_'}(male));
    p = plot(xi1,f1,xi2,f2);
    p(1).Color = '#CF0047'; p(2).Color = '#0097CD';
    p(1).LineWidth = 2; p(2).LineWidth = 2;
    legend('female','male','FontSize',12)
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
    if keep_fig == true
        title(['Distribution of ages by sex (1yr, n=',...
        num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of ages (2 yr n=', num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    figure('Name','Normplot Age - crit C')
    h1 = normplot(tbl_2yr_metabolic{:,'Age_'});
    h1(1).LineWidth = 2;h1(2).LineWidth = 2;h1(3).LineWidth = 2;
    h1(1).MarkerEdgeColor = '#E5CC00';
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age (2 year, n=', ...
        num2str(dsize),',',id_twl,'D)',fttype]);
    
    figure('Name','Normplots Age by sex')
    hold on
    hm = normplot(tbl_2yr_metabolic{:,'Age_'}(male));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_2yr_metabolic{:,'Age_'}(female));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Male','','','Female'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age by sex (2 yr, n=', ...
        num2str(dsize),',',id_twl,'D)',fttype]);
    
    % Null hypothesis: The mean strengths for the two populations are equal.
    % A p-value less than the significance level indicates that you 
    % can reject the null hypothesis. In other words, the sample provides 
    % sufficient evidence to conclude that the population means are 
    % different. Below is the output for the analysis.
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_2yr_metabolic{:,'Age_'}(male)),...
    (tbl_2yr_metabolic{:,'Age_'}(female)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_2yr_metabolic{:,'Age_'}(male)),...
    (tbl_2yr_metabolic{:,'Age_'}(female)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.C.second.tstudent(conta+1).fieldName = 'Age by sex';
    stat.C.second.tstudent(conta+1).h_ttest = httest2;
    stat.C.second.tstudent(conta+1).pval_ttest = pttest2;
    stat.C.second.tstudent(conta+1).ci_ttest = cittest2;
    stat.C.second.tstudent(conta+1).tstat = statsttest2.tstat;
    stat.C.second.tstudent(conta+1).df = statsttest2.df;
    stat.C.second.tstudent(conta+1).sd = statsttest2.sd;
    stat.C.second.tstudent(conta+1).tstudent = veredict1;
    stat.C.second.tstudent(conta+1).male = male;
    stat.C.second.tstudent(conta+1).female = female;
    stat.C.second.tstudent(conta+1).data = (tbl_2yr_metabolic(:,{'PatientCode','Sex','Age_'}));
    
    stat.C.second.Ftest(conta+1).fieldName = 'Age by sex';
    stat.C.second.Ftest(conta+1).fstat = stats_var.fstat;
    stat.C.second.Ftest(conta+1).df1 = stats_var.df1;
    stat.C.second.Ftest(conta+1).df2 = stats_var.df2;
    stat.C.second.Ftest(conta+1).h_ftest = h_var;
    stat.C.second.Ftest(conta+1).pval_ftest = p_var;
    stat.C.second.Ftest(conta+1).ci_ftest = ci_var;
    stat.C.second.Ftest(conta+1).ftest = veredict2;
    stat.C.second.Ftest(conta+1).male = male;
    stat.C.second.Ftest(conta+1).female = female;
    stat.C.second.Ftest(conta+1).data = (tbl_2yr_metabolic(:,{'PatientCode','Sex','Age_'}));
    
    tbl_cross = array2table(crosstab(female,(tbl_2yr_metabolic{:,'TWL'}))); % 1 = female
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the 2 distinct GENDERS.
    tbl_cross.Properties.RowNames = {'male','female'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association gender & TWL ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association gender & TWL';
    end
    
    stat.C.second.Fisher(conta+1).fieldName = 'Gender by class';
    stat.C.second.Fisher(conta+1).h = h;
    stat.C.second.Fisher(conta+1).pval = p;
    stat.C.second.Fisher(conta+1).OddsRatio = stats.OddsRatio;
    stat.C.second.Fisher(conta+1).ci = stats.ConfidenceInterval;
    stat.C.second.Fisher(conta+1).Fisher = veredict3;
    stat.C.second.Fisher(conta+1).contingencyTable = tbl_cross;
    stat.C.second.Fisher(conta+1).yes = 'TWL';
    stat.C.second.Fisher(conta+1).no = 'TWL';
    stat.C.second.Fisher(conta+1).data = (tbl_2yr_metabolic(:,{'PatientCode','TWL','Age_'}));
    
    clearvars male female
    
     % BMI
    success = tbl_2yr_metabolic.TWL == 1;
    fail = tbl_2yr_metabolic.TWL == 0;
    
    figure('Name','BMI 2nd yr - critD')
%     histogram(table2array(tbl_2yr_metabolic(success,{'BMI'})),...
%         'BinMethod', 'integer', 'FaceColor', '#E19B00',...
%         'EdgeColor','none')
%         hold on
%     histogram(table2array(tbl_2yr_metabolic(fail,{'BMI'})),...
%         'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
%         'EdgeColor','none');
    [f1,xi1] = ksdensity(table2array(tbl_2yr_metabolic(success,{'BMI'})));
    [f2,xi2] = ksdensity(table2array(tbl_2yr_metabolic(fail,{'BMI'})));
    p = plot(xi1,f1,xi2,f2);
    p(1).Color = '#E19B00'; p(2).Color = '#CF0047';
    p(1).LineWidth = 2; p(2).LineWidth = 2;
    legend('success','failure')
    xlabel('BMI (kg/m^2)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class among BMI min (2yr, n=', num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among BMI min (2yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    figure('Name','Normplots BMI')
    hold on
    hm = normplot(tbl_2yr_metabolic{:,'BMI'}(success));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_2yr_metabolic{:,'BMI'}(fail));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Success','','','Fail'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of BMI by class (2 yr, n=', ...
        num2str(dsize),',',id_twl,'D)',fttype]);
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_2yr_metabolic{:,'BMI'}(success)),...
    (tbl_2yr_metabolic{:,'BMI'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_2yr_metabolic{:,'BMI'}(success)),...
    (tbl_2yr_metabolic{:,'BMI'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.C.second.tstudent(conta+2).fieldName = 'BMI';
    stat.C.second.tstudent(conta+2).h_ttest = httest2;
    stat.C.second.tstudent(conta+2).pval_ttest = pttest2;
    stat.C.second.tstudent(conta+2).ci_ttest = cittest2;
    stat.C.second.tstudent(conta+2).tstat = statsttest2.tstat;
    stat.C.second.tstudent(conta+2).df = statsttest2.df;
    stat.C.second.tstudent(conta+2).sd = statsttest2.sd;
    stat.C.second.tstudent(conta+2).tstudent = veredict1;
    stat.C.second.tstudent(conta+2).sucess = success;
    stat.C.second.tstudent(conta+2).fail = fail;
    stat.C.second.tstudent(conta+2).data = (tbl_2yr_metabolic(:,{'PatientCode','TWL','BMI'}));
    
    [pWhitney,hWhitney,statsWhitney] = ranksum((tbl_2yr_metabolic{:,'BMI'}(success)),...
    (tbl_2yr_metabolic{:,'BMI'}(fail)));
    if (hWhitney == 1)
        % rejects the null hypothesis
        veredictW = 'Diff sig (median)'; 
    else
        veredictW = 'No diff sig (median)';
    end
    
    stat.C.second.Whitney(conta+1).fieldName = 'BMI';
    stat.C.second.Whitney(conta+1).h = hWhitney;
    stat.C.second.Whitney(conta+1).p = pWhitney;
    stat.C.second.Whitney(conta+1).ranksum = statsWhitney.ranksum;
    stat.C.second.Whitney(conta+1).Whitney = veredictW;
    stat.C.second.Whitney(conta+1).sucess = success;
    stat.C.second.Whitney(conta+1).fail = fail;
    stat.C.second.Whitney(conta+1).data = (tbl_2yr_metabolic(:,{'PatientCode','TWL','BMI'}));
    
    stat.C.second.Ftest(conta+2).fieldName = 'BMI';
    stat.C.second.Ftest(conta+2).fstat = stats_var.fstat;
    stat.C.second.Ftest(conta+2).df1 = stats_var.df1;
    stat.C.second.Ftest(conta+2).df2 = stats_var.df2;
    stat.C.second.Ftest(conta+2).h_ftest = h_var;
    stat.C.second.Ftest(conta+2).pval_ftest = p_var;
    stat.C.second.Ftest(conta+2).ci_ftest = ci_var;
    stat.C.second.Ftest(conta+2).ftest = veredict2;
    stat.C.second.Ftest(conta+2).sucess = success;
    stat.C.second.Ftest(conta+2).fail = fail;
    stat.C.second.Ftest(conta+2).data = (tbl_2yr_metabolic(:,{'PatientCode','TWL','BMI'}));
    
    % Age
    figure('Name','Age 2 yr - critC')
    % Uncomment in case you prefer an histogram
%     histogram(table2array(tbl_2yr_metabolic(success,{'Age_'})),...
%         'BinMethod', 'integer', 'FaceColor', '#E19B00',...
%         'EdgeColor','none')
%         hold on
%     histogram(table2array(tbl_2yr_metabolic(fail,{'Age_'})),...
%         'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
%         'EdgeColor','none');
    [f1,xi1] = ksdensity(table2array(tbl_2yr_metabolic(success,{'Age_'})));
    [f2,xi2] = ksdensity(table2array(tbl_2yr_metabolic(fail,{'Age_'})));
    p = plot(xi1,f1,xi2,f2);
    p(1).Color = '#E19B00'; p(2).Color = '#CF0047';
    p(1).LineWidth = 2; p(2).LineWidth = 2;
    legend('success','failure')
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
        hold off
        if keep_fig == true
            title(['Distribution of ages among classes (2yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of ages among classes (2yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
        end
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_2yr_metabolic{:,'Age_'}(success)),...
    (tbl_2yr_metabolic{:,'Age_'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_2yr_metabolic{:,'Age_'}(success)),...
    (tbl_2yr_metabolic{:,'Age_'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.C.second.tstudent(conta+3).fieldName = 'Age and class';
    stat.C.second.tstudent(conta+3).h_ttest = httest2;
    stat.C.second.tstudent(conta+3).pval_ttest = pttest2;
    stat.C.second.tstudent(conta+3).ci_ttest = cittest2;
    stat.C.second.tstudent(conta+3).tstat = statsttest2.tstat;
    stat.C.second.tstudent(conta+3).df = statsttest2.df;
    stat.C.second.tstudent(conta+3).sd = statsttest2.sd;
    stat.C.second.tstudent(conta+3).tstudent = veredict1;
    stat.C.second.tstudent(conta+3).sucess = success;
    stat.C.second.tstudent(conta+3).fail = fail;
    stat.C.second.tstudent(conta+3).data = (tbl_2yr_metabolic(:,{'PatientCode','TWL','BMI'}));
    
    stat.C.second.Ftest(conta+3).fieldName = 'Age and class';
    stat.C.second.Ftest(conta+3).fstat = stats_var.fstat;
    stat.C.second.Ftest(conta+3).df1 = stats_var.df1;
    stat.C.second.Ftest(conta+3).df2 = stats_var.df2;
    stat.C.second.Ftest(conta+3).h_ftest = h_var;
    stat.C.second.Ftest(conta+3).pval_ftest = p_var;
    stat.C.second.Ftest(conta+3).ci_ftest = ci_var;
    stat.C.second.Ftest(conta+3).ftest = veredict2;
    stat.C.second.Ftest(conta+3).sucess = success;
    stat.C.second.Ftest(conta+3).fail = fail;
    stat.C.second.Ftest(conta+3).data = (tbl_2yr_metabolic(:,{'PatientCode','TWL','Age_'}));
         
    clearvars success fail
    
    % Dyslipidemia
    yes = tbl_2yr_metabolic.dyslip == max(tbl_2yr_metabolic.dyslip);
    no = tbl_2yr_metabolic.dyslip == min(tbl_2yr_metabolic.dyslip);
    
    figure('Name','Dislipidemia 2 yr')
    histogram(table2array(tbl_2yr_metabolic(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_2yr_metabolic(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and dyslipedemia (2yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and dyslipedemia (2yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_2yr_metabolic{:,'TWL'}))); % 1 = dyslipidemia
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_dyslip','dyslip'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association dyslipidemia&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association dyslipidemia&TWL';
    end
   
    stat.C.second.Fisher(conta+2).fieldName = 'Dyslipidemia by class';
    stat.C.second.Fisher(conta+2).h = h;
    stat.C.second.Fisher(conta+2).pval = p;
    stat.C.second.Fisher(conta+2).OddsRatio = stats.OddsRatio;
    stat.C.second.Fisher(conta+2).ci = stats.ConfidenceInterval;
    stat.C.second.Fisher(conta+2).Fisher = veredict3;
    stat.C.second.Fisher(conta+2).yes = yes;
    stat.C.second.Fisher(conta+2).no = no;
    stat.C.second.Fisher(conta+2).contingencyTable = tbl_cross;
    stat.C.second.Fisher(conta+2).data = (tbl_2yr_metabolic(:,{'PatientCode','TWL','dyslip'}));
    
    clearvars yes no
    
    % Diabetes
    yes = tbl_2yr_metabolic.diabet == max(tbl_2yr_metabolic.diabet);
    no = tbl_2yr_metabolic.diabet == min(tbl_2yr_metabolic.diabet);
    
    figure('Name','Diabetes 2 yr')
    histogram(table2array(tbl_2yr_metabolic(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_2yr_metabolic(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
     if keep_fig == true
        title(['Distribution of success class and diabetes (2yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and diabetes (2yr, n=', num2str(dsize),',',id_twl,'D)',fttype]);
     end
     
    tbl_cross = array2table(crosstab(yes,(tbl_2yr_metabolic{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_diabetes','diabetes'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association diabetes&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association diabetes&TWL';
    end
    
    stat.C.second.Fisher(conta+3).fieldName = 'Diabetes by class';
    stat.C.second.Fisher(conta+3).h = h;
    stat.C.second.Fisher(conta+3).pval = p;
    stat.C.second.Fisher(conta+3).OddsRatio = stats.OddsRatio;
    stat.C.second.Fisher(conta+3).ci = stats.ConfidenceInterval;
    stat.C.second.Fisher(conta+3).Fisher = veredict3;
    stat.C.second.Fisher(conta+3).yes = yes;
    stat.C.second.Fisher(conta+3).no = no;
    stat.C.second.Fisher(conta+3).contingencyTable = tbl_cross;
    stat.C.second.Fisher(conta+3).data = (tbl_2yr_metabolic(:,{'PatientCode','TWL','diabet'}));
    
    clearvars yes no
    
    % Hypertension
    yes = tbl_2yr_metabolic.hypert == max(tbl_2yr_metabolic.hypert);
    no = tbl_2yr_metabolic.hypert == min(tbl_2yr_metabolic.hypert);
    
    figure('Name','Hypertension 2 yrs')
    histogram(table2array(tbl_2yr_metabolic(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_2yr_metabolic(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and hypertension (2yr, n=',num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and hypertension (2yr, n=',num2str(dsize),',',id_twl,'D)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_2yr_metabolic{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_hypert','hypert'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association hypertension&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association hypertension&TWL';
    end
    
    stat.C.second.Fisher(conta+4).fieldName = 'Hypertension by class';
    stat.C.second.Fisher(conta+4).h = h;
    stat.C.second.Fisher(conta+4).pval = p;
    stat.C.second.Fisher(conta+4).OddsRatio = stats.OddsRatio;
    stat.C.second.Fisher(conta+4).ci = stats.ConfidenceInterval;
    stat.C.second.Fisher(conta+4).Fisher = veredict3;
    stat.C.second.Fisher(conta+4).yes = yes;
    stat.C.second.Fisher(conta+4).no = no;
    stat.C.second.Fisher(conta+4).contingencyTable = tbl_cross;
    stat.C.second.Fisher(conta+4).data = (tbl_2yr_metabolic(:,{'PatientCode','TWL','hypert'}));
    
    clearvars yes no
end

if hist_dist == true
    dsize = size(tbl_2yr_bariatric,1);
    % Gender
    male = tbl_2yr_bariatric.Sex == 1;
    female = tbl_2yr_bariatric.Sex == 0;
    
    figure('Name','Gender 2yr - critD')
    histogram(tbl_2yr_bariatric{:,'TWL'}(female),'BinMethod', 'integer',...
        'FaceColor', '#CF0047','EdgeColor','none','FaceAlpha',0.4);
    hold on
    histogram(tbl_2yr_bariatric{:,'TWL'}(male),'BinMethod', 'integer',...
        'FaceColor', '#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    legend('female','male','FontSize',12)
    xlabel('Output variable (%TWL)','FontSize',16,...
        'FontWeight','bold')
    hold off
    if keep_fig == true
        title(['Distribution of TWL among gender (2yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among gender (2yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    figure('Name','Age by sex 2yr - critD')
    hold on
    histogram(tbl_2yr_bariatric{:,'Age_'}(female),'BinMethod', 'integer',...
        'FaceColor','#CF0047','FaceAlpha',0.4,'EdgeColor','none');
    histogram(tbl_2yr_bariatric{:,'Age_'}(male),'BinMethod', 'integer',...
        'FaceColor','#0097CD','FaceAlpha',0.4,'EdgeColor','none');
    hold off
    legend('female','male','FontSize',12)
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
    if keep_fig == true
        title(['Distribution of ages by sex (1yr, n=',...
        num2str(dsize),',',id_twl,')'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of ages (2 yr n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    figure('Name','Normplot Age - crit C')
    h1 = normplot(tbl_2yr_bariatric{:,'Age_'});
    h1(1).LineWidth = 2;h1(2).LineWidth = 2;h1(3).LineWidth = 2;
    h1(1).MarkerEdgeColor = '#E5CC00';
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age (1 year, n=', ...
        num2str(dsize),',',id_twl,'Dbar)',fttype]);
    
    figure('Name','Normplots Age by sex')
    hold on
    hm = normplot(tbl_2yr_bariatric{:,'Age_'}(male));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_2yr_bariatric{:,'Age_'}(female));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Male','','','Female'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of Age by sex (1 yr, n=', ...
        num2str(dsize),',',id_twl,'Dbar)',fttype]);
    
    % Null hypothesis: The mean strengths for the two populations are equal.
    % A p-value less than the significance level indicates that you 
    % can reject the null hypothesis. In other words, the sample provides 
    % sufficient evidence to conclude that the population means are 
    % different. Below is the output for the analysis.
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_2yr_bariatric{:,'Age_'}(male)),...
    (tbl_2yr_bariatric{:,'Age_'}(female)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_2yr_bariatric{:,'Age_'}(male)),...
    (tbl_2yr_bariatric{:,'Age_'}(female)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.D.second.tstudent(conta+1).fieldName = 'Age by sex';
    stat.D.second.tstudent(conta+1).h_ttest = httest2;
    stat.D.second.tstudent(conta+1).pval_ttest = pttest2;
    stat.D.second.tstudent(conta+1).ci_ttest = cittest2;
    stat.D.second.tstudent(conta+1).tstat = statsttest2.tstat;
    stat.D.second.tstudent(conta+1).df = statsttest2.df;
    stat.D.second.tstudent(conta+1).sd = statsttest2.sd;
    stat.D.second.tstudent(conta+1).tstudent = veredict1;
    stat.D.second.tstudent(conta+1).male = male;
    stat.D.second.tstudent(conta+1).female = female;
    stat.D.second.tstudent(conta+1).data = (tbl_2yr_bariatric(:,{'PatientCode','Sex','Age_'}));
    
    stat.D.second.Ftest(conta+1).fieldName = 'Age by sex';
    stat.D.second.Ftest(conta+1).fstat = stats_var.fstat;
    stat.D.second.Ftest(conta+1).df1 = stats_var.df1;
    stat.D.second.Ftest(conta+1).df2 = stats_var.df2;
    stat.D.second.Ftest(conta+1).h_ftest = h_var;
    stat.D.second.Ftest(conta+1).pval_ftest = p_var;
    stat.D.second.Ftest(conta+1).ci_ftest = ci_var;
    stat.D.second.Ftest(conta+1).ftest = veredict2;
    stat.D.second.Ftest(conta+1).male = male;
    stat.D.second.Ftest(conta+1).female = female;
    stat.D.second.Ftest(conta+1).data = (tbl_2yr_bariatric(:,{'PatientCode','Sex','Age_'}));
    
    tbl_cross = array2table(crosstab(female,(tbl_2yr_bariatric{:,'TWL'}))); % 1 = female
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the 2 distinct GENDERS.
    tbl_cross.Properties.RowNames = {'male','female'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association gender&Age ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association gender&Age';
    end
    
    stat.D.second.Fisher(conta+1).fieldName = 'Age by sex';
    stat.D.second.Fisher(conta+1).h = h;
    stat.D.second.Fisher(conta+1).pval = p;
    stat.D.second.Fisher(conta+1).OddsRatio = stats.OddsRatio;
    stat.D.second.Fisher(conta+1).ci = stats.ConfidenceInterval;
    stat.D.second.Fisher(conta+1).Fisher = veredict3;
    stat.D.second.Fisher(conta+1).contingencyTable = tbl_cross;
    stat.D.second.Fisher(conta+1).yes = 'TWL';
    stat.D.second.Fisher(conta+1).no = 'TWL';
    stat.D.second.Fisher(conta+1).data = (tbl_2yr_bariatric(:,{'PatientCode','TWL','Age_'}));
    
    clearvars male female
    
     % BMI
    success = tbl_2yr_bariatric.TWL == 1;
    fail = tbl_2yr_bariatric.TWL == 0;
    
    figure('Name','BMI 2 yr - critD')
    histogram(table2array(tbl_2yr_bariatric(success,{'BMI'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(tbl_2yr_bariatric(fail,{'BMI'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('BMI (kg/m^2)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class among BMI min (2yr, n=', num2str(dsize),',',id_twl,'Dbar)'],'FontSize',18)
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 25; %Width
        fig.Position(4) = 18.7; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class among BMI min (2yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    figure('Name','Normplots BMI')
    hold on
    hm = normplot(tbl_2yr_bariatric{:,'BMI'}(success));
    hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
    hm(1).MarkerEdgeColor = '#0097CD';
    hf = normplot(tbl_2yr_bariatric{:,'BMI'}(fail));
    hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
    hf(1).MarkerEdgeColor = '#CF0047';
    hold off
    legend({'','','Success','','','Fail'},'Location','southeast')
    xlabel('Age at surgery (1 year)','FontSize',16,'FontWeight','bold')
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Normplot of BMI by class (2 yr, n=', ...
        num2str(dsize),',',id_twl,'Dbar)',fttype]);
    
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_2yr_bariatric{:,'BMI'}(success)),...
    (tbl_2yr_bariatric{:,'BMI'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_2yr_bariatric{:,'BMI'}(success)),...
    (tbl_2yr_bariatric{:,'BMI'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.D.second.tstudent(conta+2).fieldName = 'BMI';
    stat.D.second.tstudent(conta+2).h_ttest = httest2;
    stat.D.second.tstudent(conta+2).pval_ttest = pttest2;
    stat.D.second.tstudent(conta+2).ci_ttest = cittest2;
    stat.D.second.tstudent(conta+2).tstat = statsttest2.tstat;
    stat.D.second.tstudent(conta+2).df = statsttest2.df;
    stat.D.second.tstudent(conta+2).sd = statsttest2.sd;
    stat.D.second.tstudent(conta+2).tstudent = veredict1;
    stat.D.second.tstudent(conta+2).sucess = success;
    stat.D.second.tstudent(conta+2).fail = fail;
    stat.D.second.tstudent(conta+2).data = (tbl_2yr_bariatric(:,{'PatientCode','TWL','BMI'}));
    
    [pWhitney,hWhitney,statsWhitney] = ranksum((tbl_2yr_bariatric{:,'BMI'}(success)),...
    (tbl_2yr_bariatric{:,'BMI'}(fail)));
    if (hWhitney == 1)
        % rejects the null hypothesis
        veredictW = 'Diff sig (median)'; 
    else
        veredictW = 'No diff sig (median)';
    end
    
    stat.D.second.Whitney(conta+1).fieldName = 'BMI';
    stat.D.second.Whitney(conta+1).h = hWhitney;
    stat.D.second.Whitney(conta+1).p = pWhitney;
    stat.D.second.Whitney(conta+1).ranksum = statsWhitney.ranksum;
    stat.D.second.Whitney(conta+1).Whitney = veredictW;
    stat.D.second.Whitney(conta+1).sucess = success;
    stat.D.second.Whitney(conta+1).fail = fail;
    stat.D.second.Whitney(conta+1).data = (tbl_2yr_bariatric(:,{'PatientCode','TWL','BMI'}));
    
    stat.D.second.Ftest(conta+2).fieldName = 'BMI';
    stat.D.second.Ftest(conta+2).fstat = stats_var.fstat;
    stat.D.second.Ftest(conta+2).df1 = stats_var.df1;
    stat.D.second.Ftest(conta+2).df2 = stats_var.df2;
    stat.D.second.Ftest(conta+2).h_ftest = h_var;
    stat.D.second.Ftest(conta+2).pval_ftest = p_var;
    stat.D.second.Ftest(conta+2).ci_ftest = ci_var;
    stat.D.second.Ftest(conta+2).ftest = veredict2;
    stat.D.second.Ftest(conta+2).sucess = success;
    stat.D.second.Ftest(conta+2).fail = fail;
    stat.D.second.Ftest(conta+2).data = (tbl_2yr_bariatric(:,{'PatientCode','TWL','BMI'}));
    
    % Age
    figure('Name','Age 2 yr - critD')
    histogram(table2array(tbl_2yr_bariatric(success,{'Age_'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
        hold on
    histogram(table2array(tbl_2yr_bariatric(fail,{'Age_'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('success','failure')
    xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
        hold off
        if keep_fig == true
            title(['Distribution of ages among classes (2yr, n=',...
            num2str(dsize),',',id_twl,')'],'FontSize',18)
            ax = gca; ax.FontSize = 14;
            fig = gcf;
            fig.Units = 'centimeter';
            fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
            fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
            fig.Position(3) = 20; %Width
            fig.Position(4) = 18; %Height
            saveas(gcf,[comp_foldername '/Distribution of ages among classes (2yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
        end
        
    [httest2,pttest2,cittest2,statsttest2] = ttest2((tbl_2yr_bariatric{:,'Age_'}(success)),...
    (tbl_2yr_bariatric{:,'Age_'}(fail)));
    if (httest2 == 1)
        % rejects the null hypothesis
        veredict1 = 'Diff sig (mean)'; 
    else
        veredict1 = 'No diff sig (means)';
    end
    
    [h_var,p_var,ci_var,stats_var] = vartest2((tbl_2yr_bariatric{:,'Age_'}(success)),...
    (tbl_2yr_bariatric{:,'Age_'}(fail)));
    if (h_var == 1)
        % rejects the null hypothesis
        veredict2 = 'Diff sig (var)'; 
    else
        veredict2 = 'No diff sig (var)';
    end

    stat.D.second.tstudent(conta+3).fieldName = 'Age and class';
    stat.D.second.tstudent(conta+3).h_ttest = httest2;
    stat.D.second.tstudent(conta+3).pval_ttest = pttest2;
    stat.D.second.tstudent(conta+3).ci_ttest = cittest2;
    stat.D.second.tstudent(conta+3).tstat = statsttest2.tstat;
    stat.D.second.tstudent(conta+3).df = statsttest2.df;
    stat.D.second.tstudent(conta+3).sd = statsttest2.sd;
    stat.D.second.tstudent(conta+3).tstudent = veredict1;
    stat.D.second.tstudent(conta+3).sucess = success;
    stat.D.second.tstudent(conta+3).fail = fail;
    stat.D.second.tstudent(conta+3).data = (tbl_2yr_bariatric(:,{'PatientCode','TWL','BMI'}));
    
    stat.D.second.Ftest(conta+3).fieldName = 'Age and class';
    stat.D.second.Ftest(conta+3).fstat = stats_var.fstat;
    stat.D.second.Ftest(conta+3).df1 = stats_var.df1;
    stat.D.second.Ftest(conta+3).df2 = stats_var.df2;
    stat.D.second.Ftest(conta+3).h_ftest = h_var;
    stat.D.second.Ftest(conta+3).pval_ftest = p_var;
    stat.D.second.Ftest(conta+3).ci_ftest = ci_var;
    stat.D.second.Ftest(conta+3).ftest = veredict2;
    stat.D.second.Ftest(conta+3).sucess = success;
    stat.D.second.Ftest(conta+3).fail = fail;
    stat.D.second.Ftest(conta+3).data = (tbl_2yr_bariatric(:,{'PatientCode','TWL','Age_'}));
        
    clearvars success fail
    
    % Dyslipidemia
    yes = tbl_2yr_bariatric.dyslip == max(tbl_2yr_bariatric.dyslip);
    no = tbl_2yr_bariatric.dyslip == min(tbl_2yr_bariatric.dyslip);
    
    figure('Name','Dislipidemia 2 yr')
    histogram(table2array(tbl_2yr_bariatric(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_2yr_bariatric(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and dyslipedemia (2yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and dyslipedemia (2yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_2yr_bariatric{:,'TWL'}))); % 1 = dyslipidemia
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_dyslip','dyslip'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association dyslipidemia&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association dyslipidemia&TWL';
    end
   
    stat.D.second.Fisher(conta+2).fieldName = 'Dyslipidemia by class';
    stat.D.second.Fisher(conta+2).h = h;
    stat.D.second.Fisher(conta+2).pval = p;
    stat.D.second.Fisher(conta+2).OddsRatio = stats.OddsRatio;
    stat.D.second.Fisher(conta+2).ci = stats.ConfidenceInterval;
    stat.D.second.Fisher(conta+2).Fisher = veredict3;
    stat.D.second.Fisher(conta+2).yes = yes;
    stat.D.second.Fisher(conta+2).no = no;
    stat.D.second.Fisher(conta+2).contingencyTable = tbl_cross;
    stat.D.second.Fisher(conta+2).data = (tbl_2yr_bariatric(:,{'PatientCode','TWL','dyslip'}));
    
    clearvars yes no
    
    % Diabetes
    yes = tbl_2yr_bariatric.diabet == max(tbl_2yr_bariatric.diabet);
    no = tbl_2yr_bariatric.diabet == min(tbl_2yr_bariatric.diabet);
    
    figure('Name','Diabetes 2 yr')
    histogram(table2array(tbl_2yr_bariatric(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_2yr_bariatric(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
     if keep_fig == true
        title(['Distribution of success class and diabetes (2yr, n=', num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and diabetes (2yr, n=', num2str(dsize),',',id_twl,'Dbar)',fttype]);
     end
     
    tbl_cross = array2table(crosstab(yes,(tbl_2yr_bariatric{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_diabetes','diabetes'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association diabetes&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association diabetes&TWL';
    end
    
    stat.D.second.Fisher(conta+3).fieldName = 'Diabetes by class';
    stat.D.second.Fisher(conta+3).h = h;
    stat.D.second.Fisher(conta+3).pval = p;
    stat.D.second.Fisher(conta+3).OddsRatio = stats.OddsRatio;
    stat.D.second.Fisher(conta+3).ci = stats.ConfidenceInterval;
    stat.D.second.Fisher(conta+3).Fisher = veredict3;
    stat.D.second.Fisher(conta+3).yes = yes;
    stat.D.second.Fisher(conta+3).no = no;
    stat.D.second.Fisher(conta+3).contingencyTable = tbl_cross;
    stat.D.second.Fisher(conta+3).data = (tbl_2yr_bariatric(:,{'PatientCode','TWL','diabet'}));
     
    clearvars yes no
    
    % Hypertension
    yes = tbl_2yr_bariatric.hypert == max(tbl_2yr_bariatric.hypert);
    no = tbl_2yr_bariatric.hypert == min(tbl_2yr_bariatric.hypert);
    
    figure('Name','Hypertension 2 yrs')
    histogram(table2array(tbl_2yr_bariatric(yes,{'TWL'})),...
        'BinMethod', 'integer', 'FaceColor', '#E19B00',...
        'EdgeColor','none')
    hold on
    histogram(table2array(tbl_2yr_bariatric(no,{'TWL'})),...
        'BinMethod', 'integer','FaceColor', '#CF0047','FaceAlpha',0.4,...
        'EdgeColor','none');
    legend('yes','no')
    xlabel('Success class (1=succeed,0=failed)','FontSize',16,'FontWeight','bold')
        hold off
    if keep_fig == true
        title(['Distribution of success class and hypertension (2yr, n=',num2str(dsize),',',id_twl,')'])
        ax = gca; ax.FontSize = 14;
        fig = gcf;
        fig.Units = 'centimeter';
        fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
        fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
        fig.Position(3) = 20; %Width
        fig.Position(4) = 18; %Height
        saveas(gcf,[comp_foldername '/Distribution of success class and hypertension (2yr, n=',num2str(dsize),',',id_twl,'Dbar)',fttype]);
    end
    
    tbl_cross = array2table(crosstab(yes,(tbl_2yr_bariatric{:,'TWL'}))); % 1 = diabetes
    %The rows in table correspond to the 2 distinct values in TWL, and
    %the columns correspond to the comorbidity.
    tbl_cross.Properties.RowNames = {'no_hypert','hypert'}; 
    tbl_cross.Properties.VariableNames = {'fail','success'};
    
    [h,p,stats] = fishertest(tbl_cross);
    if (h == 1)
        % Rejects null hypothesis
        veredict3 = ['Association hypertension&TWL, ',...
            tbl_cross.Properties.RowNames{2}, ' have about ',...
            num2str(stats.OddsRatio),' greater odds of being ', ...
            tbl_cross.Properties.VariableNames{2}];
    else
        veredict3 = 'No association hypertension&TWL';
    end
    
    stat.D.second.Fisher(conta+4).fieldName = 'Hypertension by class';
    stat.D.second.Fisher(conta+4).h = h;
    stat.D.second.Fisher(conta+4).pval = p;
    stat.D.second.Fisher(conta+4).OddsRatio = stats.OddsRatio;
    stat.D.second.Fisher(conta+4).ci = stats.ConfidenceInterval;
    stat.D.second.Fisher(conta+4).Fisher = veredict3;
    stat.D.second.Fisher(conta+4).yes = yes;
    stat.D.second.Fisher(conta+4).no = no;
    stat.D.second.Fisher(conta+4).contingencyTable = tbl_cross;
    stat.D.second.Fisher(conta+4).data = (tbl_2yr_bariatric(:,{'PatientCode','TWL','hypert'}));
    
    clearvars yes no
end

close all

var_names1 = tbl_in_1yr.Properties.VariableNames;
var_names2 = tbl_in_2yr.Properties.VariableNames;

%% Final Heatmaps

figure('Name','Variables values ditribution ( 1 yr - all patients)')
mat_in_1yr = table2array(tbl_in_1yr);
h6a=imagesc(mat_in_1yr(:,2:end));
colormap('bone')
% colorbar
if keep_fig == true
    ylabel('Patients (observations)','FontSize',16)
    xlabel('Features','FontSize',16)
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Features 1 yr, n=',num2str(size(tbl_in_1yr,1)),fttype]);
end

figure('Name','Variables values ditribution ( 1 yr - metabolic patients)')
mat_1yr_metabolic = table2array(tbl_1yr_metabolic);
h6b=imagesc(mat_1yr_metabolic(:,2:end));
colormap('bone')
% colorbar
if keep_fig == true
    ylabel('Patients (observations)','FontSize',16)
    xlabel('Features','FontSize',16)
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Features (metab) 1 yr, n=',num2str(size(tbl_in_1yr,1)),fttype]);
end

figure('Name','Variables values ditribution ( 1 yr - bariatric patients)')
mat_1yr_bariatric = table2array(tbl_1yr_bariatric);
h6c=imagesc(mat_1yr_bariatric(:,2:end));
colorbar
if keep_fig == true
    ylabel('Patients (observations)','FontSize',16)
    xlabel('Features','FontSize',16)
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Features (bar) 1 yr, n=',num2str(size(tbl_in_1yr,1)),fttype]);
end

figure('Name','Variables values ditribution (2 yrs)')
mat_in_2yr = table2array(tbl_in_2yr);
h7a=imagesc(mat_in_2yr(:,2:end));
colormap('bone')
% colorbar
if keep_fig == true
    ylabel('Patients (observations)','FontSize',16)
    xlabel('Features','FontSize',16)
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Features 2 yr, n=',num2str(size(tbl_in_2yr,1)),fttype]);
end

figure('Name','Variables values ditribution (2 yrs) - metab')
mat_2yr_metabolic = table2array(tbl_2yr_metabolic);
h7b=imagesc(mat_2yr_metabolic(:,2:end));
colorbar
if keep_fig == true
    ylabel('Patients (observations)','FontSize',16)
    xlabel('Features','FontSize',16)
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Features (metab) 2 yr, n=',num2str(size(tbl_in_2yr,1)),fttype]);
end

figure('Name','Variables values ditribution (2 yrs) - bar')
mat_2yr_bariatric = table2array(tbl_2yr_bariatric);
h7c=imagesc(mat_2yr_bariatric(:,2:end));
colormap('bone')
% colorbar
if keep_fig == true
    ylabel('Patients (observations)','FontSize',16)
    xlabel('Features','FontSize',16)
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
    saveas(gcf,[comp_foldername '/Features (bar) 2 yr, n=',num2str(size(tbl_in_2yr,1)),fttype]);
end

%% Correlation

R1=corrcoef(mat_in_1yr);
figure('Name','Correlation among variables 1st year')
pcolor(R1);
colorbar
if keep_fig == true
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
%     saveas(gcf,[comp_foldername '/Features 1 yr, n=',num2str(size(tbl_in_1yr,1)),fttype]);
end

R2=corrcoef(mat_in_2yr);
figure('Name','Correlation among variables 2nd year')
pcolor(R2)
colorbar
if keep_fig == true
    ax = gca; ax.FontSize = 14;
    fig = gcf;
    fig.Units = 'centimeter';
    fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
    fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
    fig.Position(3) = 20; %Width
    fig.Position(4) = 18; %Height
%     saveas(gcf,[comp_foldername '/Features 2 yr, n=',num2str(size(tbl_in_2yr,1)),fttype]);
end

%% Save final tables and matrixes

id_vars = [num2str(length(var_names1)-2),'var'];
root = 'CZE(1yr&2yr)';

if bin_class == true
    filename = sprintf('%s_%s_%s_%s',root, id_twl, id_vars);
    dat = char(datetime('now','Format','yyyy-MM-dd'));
    filename = sprintf('(%s) %s', dat, filename);
    save([comp_foldername '/' filename '.mat'],'var_names1',...
        'var_names2', 'input_1yr', 'input_2yr', ...
        'tbl_1yr_bariatric', 'tbl_1yr_metabolic', 'tbl_2yr_bariatric', ...
        'tbl_2yr_metabolic', 'tbl_in_1yr', 'tbl_in_2yr', 'var_names', ...
        'comp_foldername','Output1yr','Output2yr','ScreeningData',...
        'GeneralData','dat','stat','filename')
end
if multi_class == true
    filename = sprintf('%s_%s_%s_%s',root, 'multi',id_vars);
    dat = char(datetime('now','Format','yyyy-MM-dd'));
    filename = sprintf('(%s) %s', dat, filename);
    save([comp_foldername '/' filename '.mat'],'mat_all_1yr_norm',...
        'mat_all_2yr_norm', 'tbl_in_1yr_denorm','tbl_in_1yr_norm',...
        'tbl_in_2yr_denorm','tbl_in_2yr_norm','total_class1_1yr',...
        'total_class1_2yr','total_class2_1yr','total_class2_2yr',...
        'total_class3_1yr','total_class3_2yr','total_class4_1yr',...
        'total_class4_2yr','var_names','var_names1','var_names2',...
        'stat','filename')
end
if reg_model == true
    filename = sprintf('%s_%s_%s',root, 'reg', id_vars);
    dat = char(datetime('now','Format','yyyy-MM-dd'));
    filename = sprintf('(%s) %s', dat, filename);
    save([comp_foldername '/' filename '.mat'],'var_names1',...
        'var_names2', 'input_1yr', 'input_2yr', 'dd', ...
        'tbl_1yr_bariatric', 'tbl_1yr_metabolic', 'tbl_2yr_bariatric', ...
        'tbl_2yr_metabolic', 'tbl_in_1yr', 'tbl_in_2yr', 'var_names', ...
        'comp_foldername','stat','filename')
end

clearvars id_twl id_vars root dat h1 h h2 h3 h4 h5 h6a h6b h6c h7a h7b ...
    h7c 

%%
% _Created by Aldo Arévalo_ 