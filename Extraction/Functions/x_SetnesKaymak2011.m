%% Finding the optimal allowed percentage of missing information in a feature
% Code is based on the paper: M. Setnes and U. Kaymak. Fuzzy modeling of 
% client preference from large data sets: an application to target selection 
% in direct marketing. IEEE Transactions on Fuzzy Systems, 9(1):153–163, 
% February 2001.  
%
% Code by CF and adapted by AA
%
% The code uses InputTable as input, where each row contains the
% information about one patient and each column represents one feature.
% Every patient record has an unique identifier, in thus script called
% 'PatientNr'.

%% Determine number of NaNs per patient
function coordinates = x_SetnesKaymak2011(Input,year,comp_foldername)

InputTable=Input(:,2:end-1);

PatientNr=zeros(height(InputTable),1);
NumberOfNaNs=zeros(height(InputTable),1);
PercentageOfNaNs=zeros(height(InputTable),1);

for i=1:height(InputTable)
    PatientNr(i)=InputTable{i,1};
    NumberOfNaNs(i)=sum(isnan(InputTable{i,1:end}));
    PercentageOfNaNs(i)=(NumberOfNaNs(i)/width(InputTable))*100;
end
NaNsPerPatient=[array2table(PatientNr) array2table(NumberOfNaNs) array2table(PercentageOfNaNs)];

%% Histogram
figure(90+year)
histogram(PercentageOfNaNs,[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 ...
    75 80 85 90 95 100],'FaceColor',[0.9290 0.6940 0.1250],...
    'EdgeColor','none')
xlabel({'Percentage of NaNs'},'FontSize',16,'FontWeight','bold');
ylabel({'Number of patients'},'FontSize',16,'FontWeight','bold');
set(gca,'FontSize',18)

ax99 = gcf;
ax99.Units = 'centimeter';
ax99.Position = [15 1 28 20];

clear CaseNr NumberOfNaNs PercentageOfNaNs

print([comp_foldername,'/PatientsVSNaNs',num2str(year),'Yr'],'-deps');

ax99.PaperOrientation = 'landscape';
print([comp_foldername,'/PatientsVSNaNs',num2str(year),'Yr'],'-dpdf', '-fillpage');

%% Determine number of NaNs per variable
VariableName=InputTable.Properties.VariableNames(:,3:end)';
NumberOfNaNs=zeros(length(VariableName),1);
PercentageOfNaNs=zeros(length(VariableName),1);
[j,~]=size(VariableName);

for i=1:j
    NumberOfNaNs(i,1)=sum(isnan(InputTable{:,i+2}));
    PercentageOfNaNs(i,1)=(NumberOfNaNs(i,1)/height(InputTable))*100;
end
NaNsPerVariable=[array2table(VariableName) array2table(NumberOfNaNs) array2table(PercentageOfNaNs)];

clear j
%% Histogram
figure(100+year)
histogram(PercentageOfNaNs,[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 ...
    75 80 85 90 95 100],'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',0.7);
xlabel({'Percentage of NaNs'},'FontSize',16,'FontWeight','bold');
ylabel({'Number of variables'},'FontSize',16,'FontWeight','bold');
set(gca,'FontSize',18)

ax100 = gcf;
ax100.Units = 'centimeter';
ax100.Position = [12.6647 6.2442 31.7500 18.3444];

clear x i j VariableName NumberOfNaNs PercentageOfNaNs

print([comp_foldername,'/VariablesVSNaNs',num2str(year),'Yr'],'-deps','-tiff');

ax100.PaperOrientation = 'landscape';
print([comp_foldername,'/VariablesVSNaNs',num2str(year),'Yr'],'-dpdf', '-fillpage');

%% No. features retained
x(1)=sum(logical((NaNsPerVariable.NumberOfNaNs)==0));
x(2)=x(1)+sum(logical((NaNsPerVariable.NumberOfNaNs)==1));

for i=3:height(InputTable)+1
    x(i)=x(i-1)+sum(logical((NaNsPerVariable.NumberOfNaNs)==i));
end

%% No. of patient records retained

InputDummy=InputTable;
y=zeros(height(NaNsPerVariable),1);
NaNsPerVariable = sortrows(NaNsPerVariable,'NumberOfNaNs','ascend');

for i=1:height(NaNsPerVariable)
    d=table2array(NaNsPerVariable(i,1));
    l=logical(isnan(table2array(InputDummy(:,d)))==0);
    InputDummy=InputDummy(l,:);
    y(i,1)=height(InputDummy);
end

clear d

k=table2array(NaNsPerVariable(:,2));

z(1:(height(InputTable)+1),1)=height(InputTable);
for i=1:length(y)-1
    if k(i) ~= 0 && k(i-1) ~= 0
        z((k(i-1)+1):(k(i)+1),1)=y(i);     
    end
end

clear l i 
%% Create figure
figure(110+year)

h=0:height(InputTable);
xdimless=x/max(x);
ydimless=y/max(y);
zdimless=z/max(z);

hh = line(h,xdimless);
hh.Color = [0.8500 0.3250 0.0980];
hh.LineStyle = '--';
hh.LineWidth = 2;
xlim([0 height(InputTable)])
ylim([0 1])
hh.DisplayName = 'Features retained';
ax1=gca;
ax1.FontSize = 18;
xlabel('Number of missing data allowed in feature','FontSize',30,...
    'FontWeight','bold')

gg = line(k,ydimless);
gg.Color = [0 0.4470 0.7410];
gg.LineWidth = 2;
gg.DisplayName = 'Patient records retained';
xlim([0 height(InputTable)]);
ylim([0 1]);

% Intersection
L1 = vertcat(h,xdimless);
L2 = vertcat(k',ydimless');
coordinates = InterX(L1,L2);
x1 = [coordinates(1,1),coordinates(1,1)];
y1 = [0,coordinates(2,1)];
x2 = [0,length(h)];
y2 = [coordinates(2,1),coordinates(2,1)];
p3a = line(x1,y1,'Color',[0.4660 0.6740 0.1880],'LineStyle',':','LineWidth',2);% Vertical
p3a.DisplayName = 'Boundaries intersection';
p4a = line(x2,y2,'Color',[0.4660 0.6740 0.1880],'LineStyle',':','LineWidth',2);% Horizontal
p4a.Annotation.LegendInformation.IconDisplayStyle = 'off';

l = legend('show');
l.Location = 'best';
l.FontSize = 25;
ylabel('Fraction','FontSize',33,'FontWeight','bold')

fig = gcf;
fig.Units = 'centimeter';
fig.Position = [15 1 28 20];

print([comp_foldername,'/OptNumbPatients',num2str(year)], '-deps','-tiff')
print([comp_foldername,'/OptNumbPatientsC',num2str(year)], '-depsc','-tiff')

fig.PaperOrientation = 'landscape';
print([comp_foldername,'/OptNumbPatients',num2str(year)], '-dpdf', '-fillpage')

clear x y z xdimless ydimless zdimless h k InputDummy L1 L2
end