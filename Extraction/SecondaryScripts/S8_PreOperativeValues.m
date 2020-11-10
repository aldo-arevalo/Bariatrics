%% Code to estimate WEIGHT, BMI, HEIGHT and WAIST CIRCUMFERENCE
% Sources: DATOScreening, Physical Measures and Pre-Surgical Screening

load([subfold,'/',d,'/ImportedTables.mat'], 'PPOSdata', 'PhMsData',...
    'ScreeningData', 'total_patientsID', 'GeneralData')

A = exist('LabPackTimesPatient','var');
if A == 0
    load([subfold,'/',d,'/LabResultsTimeSeries.mat'], 'LabPackTimesPatient')
end

clearvars A

% Row's names
patients_list = total_patientsID;

% Column's names
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
        
    print([subfold,'/',d,'/Plots/Weightfinal'],'-deps')
    print([subfold,'/',d,'/Plots/Weightfinal','c'],'-depsc')
    print([subfold,'/',d,'/Plots/Weightfinal'], '-dpdf','-bestfit')
    
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
        
    print([subfold,'/',d,'/Plots/Heightfinal'],'-deps')
    print([subfold,'/',d,'/Plots/Heightfinal','c'],'-depsc')
    print([subfold,'/',d,'/Plots/Heightfinal'], '-dpdf','-bestfit')
    
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
        
    print([subfold,'/',d,'/Plots/BMIfinal'],'-deps')
    print([subfold,'/',d,'/Plots/BMIfinal','c'],'-depsc')
    print([subfold,'/',d,'/Plots/BMIfinal'], '-dpdf','-bestfit')
    
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
        
    print([subfold,'/',d,'/Plots/WCfinal'],'-deps')
    print([subfold,'/',d,'/Plots/WCfinal','c'],'-depsc')
    print([subfold,'/',d,'/Plots/WCfinal'], '-dpdf','-bestfit')
    
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
                        weight_max.Max(i) = max(a);
%                         weight_max.Max(i) = min(a);
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

%% Standarize BRI
 % Gather all BRI values
BRIall=vertcat(PreOpMedian.BRI,PreOpMean.BRI,PreOpMin.BRI,PreOpMax.BRI);
BRIall(isnan(BRIall))=[]; BRIall(~any(BRIall,2))=[];
[D,PD] = allfitdist(BRIall);
    
mBRI = NaN(3,1);
sBRI = NaN(3,1);

DBRI = D; PDBRI = PD;
    
for i=1:3
    s = std(PD{i}); sBRI(i,1) = s;
    m = mean(PD{i}); mBRI(i,1) = m;
    figure('Name',['BRI ',D(i).DistName])
    nbins = max(min(length(BRIall)./10,100),50);
    xi = linspace(min(BRIall),max(BRIall),nbins);
    dx = mean(diff(xi));
    xi2 = linspace(min(BRIall),max(BRIall),nbins*10)';
    fi = histc(BRIall,xi-dx);
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
    xlabel('BRI','FontSize',18,'FontWeight','bold') % x-axis label
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

        print([subfold,'/',d,'/Plots/BRI',D(i).DistName],'-deps')
        print([subfold,'/',d,'/Plots/BRI',D(i).DistName,'c'],'-depsc')
%         print([subfold,'/',d,'/Plots/BRI',D(i).DistName], '-dpdf','-bestfit')
end
    
close all
clearvars ax b1 D dx fi fig h1 h2 h3 h4 h5 i m nbins p1 PD s xi xi2 xL ys

% Standarize values
s = ((mBRI(1,1)+2*sBRI(1,1)) - (mBRI(1,1)-2*sBRI(1,1)))/4; % sigma = (RefUp - RefLow)/4
m = mBRI(1,1);        % mu = RefLow + 2s
for pp=1:size(PreOpMedian,1)
    PreOpMedian(pp,{'BRI'}) = array2table((table2array(PreOpMedian(pp,{'BRI'}))-m)/s);
    PreOpMean(pp,{'BRI'}) = array2table((table2array(PreOpMean(pp,{'BRI'}))-m)/s);
    %PreOpMin(pp,{'BRI'}) = array2table((table2array(PreOpMin(pp,{'BRI'}))-m)/s);
    %PreOpMax(pp,{'BRI'}) = array2table((table2array(PreOpMax(pp,{'BRI'}))-m)/s);
end
clearvars s m pp

%% Standarize ABSI
 % Gather all ABSI values
ABSIall=vertcat(PreOpMedian.ABSI,PreOpMean.ABSI,PreOpMin.ABSI,PreOpMax.ABSI);
ABSIall(isnan(ABSIall))=[]; ABSIall(~any(ABSIall,2))=[];
[D,PD] = allfitdist(ABSIall);
    
mABSI = NaN(3,1);
sABSI = NaN(3,1);

DABSI = D; PDABSI = PD;
    
for i=1:3
    s = std(PD{i}); sABSI(i,1) = s;
    m = mean(PD{i}); mABSI(i,1) = m;
    figure('Name',['ABSI ',D(i).DistName])
    nbins = max(min(length(ABSIall)./10,100),50);
    xi = linspace(min(ABSIall),max(ABSIall),nbins);
    dx = mean(diff(xi));
    xi2 = linspace(min(ABSIall),max(ABSIall),nbins*10)';
    fi = histc(ABSIall,xi-dx);
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
    xlabel('ABSI','FontSize',18,'FontWeight','bold') % x-axis label
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

%         print([subfold,'/',d,'/Plots/ABSI',D(i).DistName],'-deps')
%         print([subfold,'/',d,'/Plots/ABSI',D(i).DistName,'c'],'-depsc')
%         print([subfold,'/',d,'/Plots/ABSI',D(i).DistName], '-dpdf','-bestfit')
end
    
close all
clearvars ax b1 D dx fi fig h1 h2 h3 h4 h5 i m nbins p1 PD s xi xi2 xL ys

% Standarize values
s = ((mABSI(1,1)+2*sABSI(1,1)) - (mABSI(1,1)-2*sABSI(1,1)))/4; % sigma = (RefUp - RefLow)/4
m = mABSI(1,1);        % mu = RefLow + 2s
for pp=1:size(PreOpMedian,1)
    PreOpMedian(pp,{'ABSI'}) = array2table((table2array(PreOpMedian(pp,{'ABSI'}))-m)/s);
    PreOpMean(pp,{'ABSI'}) = array2table((table2array(PreOpMean(pp,{'ABSI'}))-m)/s);
    %PreOpMin(pp,{'ABSI'}) = array2table((table2array(PreOpMin(pp,{'ABSI'}))-m)/s);
    %PreOpMax(pp,{'ABSI'}) = array2table((table2array(PreOpMax(pp,{'ABSI'}))-m)/s);
end
clearvars s m pp

%% Standarize TBFM
 % Gather all TBFM values
% TBFMall=vertcat(PreOpMedian.TBFM,PreOpMean.TBFM,PreOpMin.TBFM,PreOpMax.TBFM);
TBFMall=PreOpLast.TBFM;
TBFMall(isnan(TBFMall))=[]; TBFMall(~any(TBFMall,2))=[];
[D,PD] = allfitdist(TBFMall);
    
mTBFM = NaN(3,1);
sTBFM = NaN(3,1);

DTBFM = D; PDTBFM = PD;
    
for i=1:3
    s = std(PD{i}); sTBFM(i,1) = s;
    m = mean(PD{i}); mTBFM(i,1) = m;
    figure('Name',['TBFM ',D(i).DistName])
    nbins = max(min(length(TBFMall)./10,100),50);
    xi = linspace(min(TBFMall),max(TBFMall),nbins);
    dx = mean(diff(xi));
    xi2 = linspace(min(TBFMall),max(TBFMall),nbins*10)';
    fi = histc(TBFMall,xi-dx);
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
    xlabel('TBFM','FontSize',18,'FontWeight','bold') % x-axis label
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

%         print([subfold,'/',d,'/Plots/TBFM',D(i).DistName],'-deps')
%         print([subfold,'/',d,'/Plots/TBFM',D(i).DistName,'c'],'-depsc')
%         print([subfold,'/',d,'/Plots/TBFM',D(i).DistName], '-dpdf','-bestfit')
end
    
close all
clearvars ax b1 D dx fi fig h1 h2 h3 h4 h5 i m nbins p1 PD s xi xi2 xL ys

% Standarize values
s = ((mTBFM(1,1)+2*sTBFM(1,1)) - (mTBFM(1,1)-2*sTBFM(1,1)))/4; % sigma = (RefUp - RefLow)/4
m = mTBFM(1,1);        % mu = RefLow + 2s
for pp=1:size(PreOpMedian,1)
    PreOpMedian(pp,{'TBFM'}) = array2table((table2array(PreOpMedian(pp,{'TBFM'}))-m)/s);
    PreOpMean(pp,{'TBFM'}) = array2table((table2array(PreOpMean(pp,{'TBFM'}))-m)/s);
    %PreOpMin(pp,{'TBFM'}) = array2table((table2array(PreOpMin(pp,{'TBFM'}))-m)/s);
    %PreOpMax(pp,{'TBFM'}) = array2table((table2array(PreOpMax(pp,{'TBFM'}))-m)/s);
    PreOpLast(pp,{'TBFM'}) = array2table((table2array(PreOpLast(pp,{'TBFM'}))-m)/s);
end
clearvars s m pp

%% Bar plot
% Frequency of variable measurements

% Find how many questionIDs have each patient
f=sum(~ismissing(PreOpNM{:,2:end},0));

var_gradient = ["#E5CC00", "#D5AA1F", "#C6883F", "#B7665F", "#A7447F",...
    "#98219F", "#8900BF"];

var_gradient_rgb = zeros(length(var_gradient),3);
    
for i = 1:length(var_gradient)
    color_str = var_gradient{i};
    color_bar = sscanf(color_str(2:end),'%2x%2x%2x',[1 3])/255;
    var_gradient_rgb(i,:) = color_bar;
end

clearvars i color_bar color_str

figure('Name','Frequency of each measurement per patient');
g = bar(f,'EdgeColor','none','FaceColor','flat','FaceAlpha',0.8);
g.CData = var_gradient_rgb;
hline = refline([0 1000]);
hline.Color = 'w'; hline.LineWidth = 4;
xlabel(['Measurements, n=',num2str(length(f))],'FontSize',16,...
    'FontWeight','bold');
ylabel('Frequency','FontSize',16,'FontWeight','bold');
xticklabels(PreOpNM.Properties.VariableNames(2:end))
fig = gcf; ax = gca;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 28; %Width
fig.Position(4) = 20; %Height
ax.FontSize = 14;
print([subfold,'/',d,'/Plots/FreqPreOpVal'],'-deps')
print([subfold,'/',d,'/Plots/FreqPreOpValc'],'-depsc')
fig.PaperOrientation = 'landscape';
print([subfold,'/',d,'/Plots/FreqPreOpVal'], '-dpdf','-bestfit')

clearvars hline

%% Histogram per variable
n_vars = length(f);
var_list ={'A) Weight','B) BMI','C) Height','D) WC','E) BRI','F) ABSI','G) TBFM'};
units = {'kg', 'kg/m^2', 'cm', 'cm', 'Adimensional', 'Adimensional', 'kg'};
subplot_idx = 0;

figure('Name','Indexes and physical measures (max vs. min)');
for p = 1:n_vars
    subplot_idx = subplot_idx + 1;
    subplot(2,4,subplot_idx)
    if subplot_idx == 6
        q(p) = histogram(PreOpMax{:, p+1},10);
        q(p).EdgeColor = 'none';
        q(p).FaceColor = '#E5CD00';
        q(p).FaceAlpha = 0.7;
        hold on
        h(p) = histogram(PreOpMin{:, p+1},10);
        h(p).EdgeColor = 'none';
        h(p).FaceColor = '#D20028';
        h(p).FaceAlpha = 0.4;
        % xlim([0 Inf])
        xlabel(cellstr(units(p)),'FontSize',16,'FontWeight','bold');
        title(cellstr(var_list(p)),'FontSize',18,'FontWeight','bold');
        legend('Maximum','Minimum');
        hold off
    else
        q(p) = histogram(PreOpMax{:, p+1},'BinMethod','integer');
        q(p).EdgeColor = 'none';
        q(p).FaceColor = '#E5CD00';
        q(p).FaceAlpha = 0.7;
        hold on
        h(p) = histogram(PreOpMin{:, p+1},'BinMethod','integer');
            %'FaceAlpha',0.2);
        h(p).EdgeColor = 'none';
        h(p).FaceColor = '#D20028';
        h(p).FaceAlpha = 0.4;
        % xlim([0 Inf])
        xlabel(cellstr(units(p)),'FontSize',16,'FontWeight','bold');
        title(cellstr(var_list(p)),'FontSize',18,'FontWeight','bold');
        legend('Maximum','Minimum');
        hold off
    end
end

fig = gcf; ax = gca;
fig.Units = 'centimeter';
fig.Position(1) = 15; %Distance from the left edge of the primary display to the inner left edge of the figure window.
fig.Position(2) = 1; %Distance from the bottom edge of the primary display to the inner bottom edge of the figure window. 
fig.Position(3) = 28; %Width
fig.Position(4) = 20; %Height
ax.FontSize = 14;

print([subfold,'/',d,'/Plots/MinMaxPreOpVal'],'-deps')
print([subfold,'/',d,'/Plots/MinMaxPreOpValc'],'-depsc')
fig.PaperOrientation = 'landscape';
print([subfold,'/',d,'/Plots/MinMaxPreOpVal'], '-dpdf','-bestfit')

clearvars n_vars f p var_list units subplot_idx;
close all

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