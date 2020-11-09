function plotmfs_Zscore(mfs,opt,dom,dom2,features,c_display)
% PLOTMFS    Plot membership functions.
%    PLOTMFS(MFS,OPT,FIG) plots membership functions. MFS can be a numerical
%    array, a cell array or a fuzzy model structure FM. The OPT parameter
%    can be used to specify the PLOT options (optional). The FIG parameter 
%    (optional) specifies the figure where to begin  with the plots. 
%    The default is 1. 

% Copyright (c) Robert Babuska, 1998.

% dom = dataset original (double)
% dom2 = dataset normalized (double)

FontSize = 14;
if nargin < 2, opt = ''; elseif isempty(opt), opt = ''; end;
if nargin < 3, dom = []; fig = 1; 
elseif length(dom) == 1, fig = dom; dom = []; 
else fig = 1; end;

cgreen = [0.4660, 0.6740, 0.1880];
cblue = [0.3010, 0.7450, 0.9330];
cyellow = [0.9290, 0.6940, 0.1250];

if isnumeric(mfs) && ~isempty(mfs)
   if isempty(dom) 
     dom = [min(mfs(:,2)) max(mfs(:,5))];
     dom = (dom(1) : (dom(2)-dom(1))/500 : dom(2))';
   end;
   if strcmp(features,'Hematos')
       x = [min(dom2) -1.6 -1.6 min(dom2)]; y = [0 0 1.05 1.05]; %hematos low
       patch(x,y,cblue,'FaceAlpha',.2)
       x = [max(dom2) 2.3 2.3 max(dom2)]; y = [0 0 1.05 1.05]; %hematos above
       patch(x,y,cyellow,'FaceAlpha',.2)
       x = [-1.6 2.3 2.3 -1.6]; y = [0 0 1.05 1.05]; %hematos normal
       patch(x,y,cgreen,'FaceAlpha',.2)
   end
   if strcmp(features,'p[Bili]')
       x = [min(dom2) -1.29 -1.29 min(dom2)]; y = [0 0 1.05 1.05];%Bili low
       patch(x,y,cblue,'FaceAlpha',.2)
       x = [-1.29 max(dom2) max(dom2) -1.29]; y = [0 0 1.05 1.05];%Bili normal
       patch(x,y,cgreen,'FaceAlpha',.2)
   end
   if strcmp(features,'p[Ch]')
        x = [min(dom2) 2.152 2.152 min(dom2)]; y = [0 0 1.05 1.05];%normal chol
        patch(x,y,cgreen,'FaceAlpha',.2)
        x = [max(dom2) 2.152 2.152 max(dom2)]; y = [0 0 1.05 1.05];%high chol
        patch(x,y,cyellow,'FaceAlpha',.2)
   end
   if strcmp(features,'BRImax')
       x = [min(dom2) -1.26 -1.26 min(dom2)]; y = [0 0 1.05 1.05];%Obsese
       patch(x,y,cblue,'FaceAlpha',.2)
       x = [max(dom2) -1.26 -1.26 max(dom2)]; y = [0 0 1.05 1.05];
       % Extreme obese
       patch(x,y,cyellow,'FaceAlpha',.2)
   end
   for q=1:size(mfs,1)
       hold on
       %dom2 = original dataset
       %plot(sort(dom),mgrade(sort(dom2),mfs(q,:)),opt{q,1},'LineWidth',2,'DisplayName',c_display{q,1});
       plot(sort(dom2),mgrade(sort(dom2),mfs(q,:)),opt{q,1},'LineWidth',2,'DisplayName',c_display{q,1});
   end
   hold off
   %set(gca,'YLim',[0 1.05],'XLim',[min(dom) max(dom)],'FontSize',FontSize,'FontWeight','bold');
   set(gca,'YLim',[0 1.05],'XLim',[min(dom2) max(dom2)],'FontSize',FontSize,'FontWeight','bold');
   %set(gca,'YLim',[0 1.05],'XLim',[0 1],'FontSize',FontSize,'FontWeight','bold');
   xlabel(regexprep(features,'_',' '),'FontSize',FontSize,'FontWeight','bold');
   ylabel('\mu','FontSize',16);
   %grid on
   
   %l = legend('show');
   %l.Location = 'west';
   %l.Location = 'best';
   %l.Orientation = 'horizontal';
   %l.FontSize = 14;
   
elseif isstruct(mfs)
   FM = mfs;
	vars = antename(FM);
   for i = 1 : FM.no
      figure(fig+i-1); clf; set(gcf,'numbertitle','off','name',['Membership functions (' FM.OutputName{i} ')']);
      if mfs.ante(i) == 2
         [~,nvar] = size(FM.rls{i});% size(FM.mfs{i});
         antepos = 1;
         antediffers = getregresdiff(FM,'ante',i);
         for k = 1 : nvar
            switch nvar
               case 1, subplot(1,1,1);
               case {2,3}, subplot(nvar,1,k);
               case {4,5,6,7,8}, subplot(ceil(nvar/2),2,k);
               otherwise, subplot(ceil(nvar/4),4,k)
            end;      
            
            while ~antediffers(1,antepos) && antepos < length(antediffers) 
               antepos = antepos + 1;
            end   
   
            dom = [FM.zmin{i}(antepos) FM.zmax{i}(antepos)];
            dom = (dom(1) : (dom(2)-dom(1))/500 : dom(2))';
            antepos = antepos + 1;
            
            plotmfs(FM.mfs{i}{k},opt,dom);
            set(gca,'FontSize',FontSize)
            xlabel(vars{i}{k},'FontSize',9);
            ylabel('\mu','FontSize',9);
            set(gca,'xlim',[FM.zmin{i}(antepos-1) FM.zmax{i}(antepos-1)]);
         end;   
      else
         [c,nvar] = size(FM.V{i});
         if c > 1
            for k = 1 : nvar            
               dom = [FM.zmin{i}(k) FM.zmax{i}(k)];
               dom = (dom(1) : (dom(2)-dom(1))/500 : dom(2))';
               for j = 1 : c
                  xx = ones(501,1)*FM.V{i}(j,:);
                  xx(:,k) = dom;
                  ff = fgrade(FM,xx,i);
                  f(:,j) = ff(:,j);
               end;
               switch nvar
               case 1, subplot(1,1,1);
               case {2,3}, subplot(nvar,1,k);
               otherwise, subplot(ceil(nvar/2),2,k);
               end;      
               plot(dom,f,opt);
               set(gca,'YLim',[0 1]);
               set(gca,'FontSize',FontSize)
               xlabel(vars{i}{k},'FontSize',9);
               ylabel('\mu','FontSize',9);
            end;
         else
            text(.5,.5,'no membership functions available ...','horizontalalignment','center');
         end;
      end;   
   end;   
elseif iscell(mfs)
   n = length(mfs);
   for i = 1 : n
      switch n
         case 1, subplot(1,1,1);
         case {2,3}, subplot(n,1,i);
         case {4,5,6,7,8}, subplot(ceil(n/2),2,i);
         otherwise, subplot(ceil(n/4),4,i)
         %otherwise, subplot(ceil(n/2),2,i)
      end;      
      plotmfs_Zscore(mfs{i},opt,dom(:,i),dom2(:,i),features{i},c_display);
   end;
end;

