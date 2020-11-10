%% General data

% load data
load([subfold,'/',d,'/RawData/GeneralDataRAW.mat'])

comp_foldername = [subfold,'/',d,'/Plots'];

dsize = size(GeneralData,1); fttype = '.png';

% Age
male = GeneralData.Sex == 1;
female = GeneralData.Sex == 0;

% Distributions
% If h = 1, this indicates the rejection of the null hypothesis at the ...
% alpha significance level. It has normal distribution
% The alternative hypothesis assumes that some difference exists between
% the true mean (?) and the comparison value (m0), whereas the null 
% hypothesis assumes that no difference exists.
[h_all,pval_all,jbstat_all,critval_all] = jbtest(GeneralData.Age);
[h_m,pval_m,jbstat_m,critval_m] = jbtest(GeneralData{:,'Age'}(male));
[h_f,pval_f,jbstat_f,critval_f] = jbtest(GeneralData{:,'Age'}(female));

figure('Name','Normplot Age')
h1 = normplot(GeneralData.Age);
h1(1).LineWidth = 2;h1(2).LineWidth = 2;h1(3).LineWidth = 2;
h1(1).MarkerEdgeColor = '#E5CC00';
xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
ax = gca; ax.FontSize = 14;
fig = gcf;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 20; %Width
fig.Position(4) = 18; %Height
saveas(gcf,[comp_foldername '/Normplot of age (raw) (n=', num2str(dsize),')',fttype]);


figure('Name','Normplots Age by sex')
hold on
hm = normplot(GeneralData{:,'Age'}(male));
hm(1).LineWidth = 2;hm(2).LineWidth = 2;hm(3).LineWidth = 2;
hm(1).MarkerEdgeColor = '#0097CD';
hf = normplot(GeneralData{:,'Age'}(female));
hf(1).LineWidth = 2;hf(2).LineWidth = 2;hf(3).LineWidth = 2;
hf(1).MarkerEdgeColor = '#CF0047';
hold off
legend({'','','Male','','','Female'},'Location','southeast')
xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
ax = gca; ax.FontSize = 14;
fig = gcf;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 20; %Width
fig.Position(4) = 18; %Height
saveas(gcf,[comp_foldername '/Normplot of ages by sex (raw) (n=', num2str(dsize),')',fttype]);


[hvartest2,pvartest2,civartest2,statsvartest2] = vartest2((GeneralData{:,'Age'}(male)),...
    (GeneralData{:,'Age'}(female)));
[httest2,pttest2,cittest2,statsttest2] = ttest2((GeneralData{:,'Age'}(male)),...
    (GeneralData{:,'Age'}(female)));

figure('Name','Age Raw by sex')
hold on
histogram(GeneralData{:,'Age'}(female),'BinMethod', 'integer',...
    'FaceColor','#CF0047','FaceAlpha',0.4,'EdgeColor','none');
histogram(GeneralData{:,'Age'}(male),'BinMethod', 'integer',...
    'FaceColor','#0097CD','FaceAlpha',0.4,'EdgeColor','none');
hold off
legend('female','male','FontSize',12)
xlabel('Age at surgery (years)','FontSize',16,'FontWeight','bold')
title(['Distribution of ages by sex (Raw, n=',...
num2str(dsize),')'],'FontSize',18)
ax = gca; ax.FontSize = 14;
fig = gcf;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 20; %Width
fig.Position(4) = 18; %Height
saveas(gcf,[comp_foldername '/Distribution of ages (raw) (n=', num2str(dsize),')',fttype]);

figure('Name','Age by sex (raw)')
[f1,xi1] = ksdensity(GeneralData{:,'Age'}(female));
[f2,xi2] = ksdensity(GeneralData{:,'Age'}(male));
p = plot(xi1,f1,xi2,f2);
p(1).Color = '#CF0047'; p(2).Color = '#0097CD';
p(1).LineWidth = 2; p(2).LineWidth = 2;
legend('female','male','FontSize',12)
xlabel('Age at surgery (raw)','FontSize',16,'FontWeight','bold')
title(['Distribution of ages by sex (raw, n=',...
num2str(dsize),')'],'FontSize',18)
ax = gca; ax.FontSize = 14;
fig = gcf;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 20; %Width
fig.Position(4) = 18; %Height
saveas(gcf,[comp_foldername '/Distribution of ages (Raw n=', num2str(dsize),'D)',fttype]);

clearvars female male ax fig p xi1 xi2 hf hm h1 f1 f2

close all

%% Comorbidities

% Load data
load([subfold,'/',d,'/RawData/ScreeningDataRAW.mat'])

% Diabetes
diabet = ScreeningData.diabet == 1;
nodiabet = ScreeningData.diabet == 0;

% Dyslipidemia
dyslip = ScreeningData.dyslip == 1;
nodyslip = ScreeningData.dyslip == 0;

% Hypertension
hyper = ScreeningData.hypert == 1;
nohyper = ScreeningData.hypert == 0;

%% Hypothesis testing final datasets
load([subfold,'/',d,'/Input_tables_(25%TWL)/',filename,'.mat'],'stat')
load([subfold,'/',d,'/Classifiers',char(string(year(datetime(date)))),...
    '+-3monthsRAW.mat'], 'Classifiers1yr','Classifiers2yr')

%% Hypothesis testing final datasets 1st year
both_1yr = table2array(stat.B.first.Fisher(1).data(:,'TWL'));
metabolic_1yr = table2array(stat.C.first.Fisher(1).data(:,'TWL'));
bariatric_1yr = table2array(stat.D.first.Fisher(1).data(:,'TWL'));
tbl_pval_1yr = zeros(4,4); tbl_pval_1yr = array2table(tbl_pval_1yr);
tbl_pval_1yr.Properties.RowNames = {'Both','Bariatric','Metabolic',...
    'RAW_1yr'}; 
tbl_pval_1yr.Properties.VariableNames = {'Both','Bariatric','Metabolic',...
    'RAW_1yr'};
tbl_h_1yr = tbl_pval_1yr;

% Create for the 2nd year
tbl_pval_2yr = tbl_pval_1yr; tbl_h_2yr = tbl_pval_1yr;

% Estimate TWL for Raw data
RAW_1yr = innerjoin(GeneralData,Classifiers1yr(:,{'PatientCode','TWL'}),...
    'LeftKeys','PatientNr','RightKeys','PatientCode');
RAW_1yr(:,'twl') = array2table(table2array(RAW_1yr(:,'TWL')) >= 0.25);
RAW_1yr(isnan(table2array(RAW_1yr(:,'TWL'))),:)=[];

% Rows: datasets, Columns: Success (class 1), Fail (class 0)
tbl_cross = zeros(4,2);
tbl_cross(1,1) = sum(both_1yr); 
tbl_cross(1,2) = length(both_1yr)-sum(both_1yr);
tbl_cross(2,1) = sum(metabolic_1yr);
tbl_cross(2,2) = length(metabolic_1yr)-sum(metabolic_1yr);
tbl_cross(3,1) = sum(bariatric_1yr);
tbl_cross(3,2) = length(bariatric_1yr)-sum(bariatric_1yr);
tbl_cross(4,1) = sum(single(table2array(RAW_1yr(:,'twl'))));
tbl_cross(4,2) = size(RAW_1yr,1) - sum(single(table2array(RAW_1yr(:,'twl'))));
tbl_cross1yr = array2table(tbl_cross);
tbl_cross1yr.Properties.RowNames = {'Both','Bariatric','Metabolic','RAW_1yr'}; 
tbl_cross1yr.Properties.VariableNames = {'Success','Fail'};

[h,p,~] = fishertest(tbl_cross1yr({'Both','RAW_1yr'},:));
tbl_pval_1yr('Both','RAW_1yr') = array2table(p); tbl_h_1yr('Both','RAW_1yr') = array2table(h);

[h,p,~] = fishertest(tbl_cross1yr({'Metabolic','RAW_1yr'},:));
tbl_pval_1yr('Metabolic','RAW_1yr') = array2table(p); tbl_h_1yr('Metabolic','RAW_1yr') = array2table(h);

[h,p,~] = fishertest(tbl_cross1yr({'Metabolic','Both'},:));
tbl_pval_1yr('Metabolic','Both') = array2table(p); tbl_h_1yr('Metabolic','Both') = array2table(h);

[h,p,~] = fishertest(tbl_cross1yr({'Bariatric','RAW_1yr'},:));
tbl_pval_1yr('Bariatric','RAW_1yr') = array2table(p); tbl_h_1yr('Bariatric','RAW_1yr') = array2table(h);

[h,p,~] = fishertest(tbl_cross1yr({'Bariatric','Both'},:));
tbl_pval_1yr('Bariatric','Both') = array2table(p); tbl_h_1yr('Bariatric','Both') = array2table(h);

[h,p,~] = fishertest(tbl_cross1yr({'Bariatric','Metabolic'},:));
tbl_pval_1yr('Bariatric','Metabolic') = array2table(p); tbl_h_1yr('Bariatric','Metabolic') = array2table(h);


%% Hypothesis testing final datasets 2nd year
both_2yr = table2array(stat.B.second.Fisher(1).data(:,'TWL'));
metabolic_2yr = table2array(stat.C.second.Fisher(1).data(:,'TWL'));
bariatric_2yr = table2array(stat.D.second.Fisher(1).data(:,'TWL'));

% Estimate TWL for Raw data
RAW_2yr = innerjoin(GeneralData,Classifiers2yr(:,{'PatientCode','TWL'}),...
    'LeftKeys','PatientNr','RightKeys','PatientCode');
RAW_2yr(:,'twl') = array2table(table2array(RAW_2yr(:,'TWL')) >= 0.25);
RAW_2yr(isnan(table2array(RAW_2yr(:,'TWL'))),:)=[];

% Rows: datasets, Columns: Success (class 1), Fail (class 0)
tbl_cross = zeros(4,2);
tbl_cross(1,1) = sum(both_2yr); 
tbl_cross(1,2) = length(both_2yr)-sum(both_2yr);
tbl_cross(2,1) = sum(metabolic_2yr);
tbl_cross(2,2) = length(metabolic_2yr)-sum(metabolic_2yr);
tbl_cross(3,1) = sum(bariatric_2yr);
tbl_cross(3,2) = length(bariatric_2yr)-sum(bariatric_2yr);
tbl_cross(4,1) = sum(single(table2array(RAW_2yr(:,'twl'))));
tbl_cross(4,2) = size(RAW_2yr,1) - sum(single(table2array(RAW_2yr(:,'twl'))));
tbl_cross2yr = array2table(tbl_cross);
tbl_cross2yr.Properties.RowNames = {'Both','Bariatric','Metabolic','RAW_2yr'}; 
tbl_cross2yr.Properties.VariableNames = {'Success','Fail'};

[h,p,~] = fishertest(tbl_cross2yr({'Both','RAW_2yr'},:));
tbl_pval_2yr('Both','RAW_2yr') = array2table(p); tbl_h_2yr('Both','RAW_2yr') = array2table(h);

[h,p,~] = fishertest(tbl_cross2yr({'Metabolic','RAW_2yr'},:));
tbl_pval_2yr('Metabolic','RAW_2yr') = array2table(p); tbl_h_2yr('Metabolic','RAW_2yr') = array2table(h);

[h,p,~] = fishertest(tbl_cross2yr({'Metabolic','Both'},:));
tbl_pval_2yr('Metabolic','Both') = array2table(p); tbl_h_2yr('Metabolic','Both') = array2table(h);

[h,p,~] = fishertest(tbl_cross2yr({'Bariatric','RAW_2yr'},:));
tbl_pval_2yr('Bariatric','RAW_2yr') = array2table(p); tbl_h_2yr('Bariatric','RAW_2yr') = array2table(h);

[h,p,~] = fishertest(tbl_cross2yr({'Bariatric','Both'},:));
tbl_pval_2yr('Bariatric','Both') = array2table(p); tbl_h_2yr('Bariatric','Both') = array2table(h);

[h,p,~] = fishertest(tbl_cross2yr({'Bariatric','Metabolic'},:));
tbl_pval_2yr('Bariatric','Metabolic') = array2table(p); tbl_h_2yr('Bariatric','Metabolic') = array2table(h);

%%
% _Created by Aldo Arévalo_ 