%% Reads in two flat files and generates the boxplot comparing
% simulation data against the GMPEs. Mostly unmodified version
% of the script from Ronnie Kamai

ymin = 0.0003;
ymax = 20000*ymin;

% Name = ['-gmpe' num2str(Mag*10) num2str(Dist) ];
Sims = csvread(SimsFile);
GMPEnum = csvread(GMPEFile);
GMPELabels = textread(GMPELabels, '%s', 'whitespace', ',');

T = [0.010,0.020,0.029,0.04,0.05,0.06,0.075,0.1,0.15,0.2,0.26,0.3,0.4,0.5,...
      0.6,0.75,1,1.5,2,3,4,5,6,7.5,10];

% Find column numbers corresponding to T
ind1 = zeros(1,size(GMPEnum,2));
ind2 = zeros(1,size(Sims,2));
for j=1:length(T);
    Tt = T(j);
    indT1 = GMPEnum(1,:)==Tt;
    indT2 = Sims(1,:)==Tt;
    ind1 = ind1 + indT1;
    ind2 = ind2 + indT2;
end

colsGMPET = logical(ind1);
colsSimT = logical(ind2);

%% Build the GMPE matrices
GMPET = GMPEnum(1,colsGMPET);
GMPEindex = GMPEnum(:,1);
for i=1:NumGMPE;
  GMPEdata(:,:,i) = GMPEnum(GMPEindex==i,colsGMPET);
end

%% Build the Simulation results matrices: 
% X is Period, Y is RotD050. get mean and Std.
X = Sims(1,colsSimT);
Y = Sims(2:end,colsSimT);
meanY = exp(mean(log(Y)));
stdY = std(Y)/sqrt(length(Y));

%% Colors
RGB_xls=[[0.96863  0.58824  0.27451];...  %9 line colors in 'excel' colormap
    [0.50196  0.39216  0.77647];...
    [0.60784  0.73333  0.34902];...
    [0.75294  0.31373  0.30196];...
    [0.58039  0.54118  0.32941];...
    [0.29412  0.67451  0.77647];...
    [0.30980  0.50588  0.74118];...
    [0.12157  0.28627  0.49020];...
    [0.45  0.45  0.45]];

%% Markers
PLOT_markers=['o', 's', 'd', 'x', '*', '+', 'v', 'h', 'p'];

%% Set up plot info

meansim = meanY;
stdsim = stdY;
Ysim = Y;
boxposition = X';
boxwidth = sortrows(1.5*10.^(log10(X')-1)); %for 'traditional' plotstyle

for i=1:NumGMPE;
YGMPE(:,i) = exp(mean(log(GMPEdata(:,:,i))));
end

YGMPE = YGMPE';
GMPEAve = exp(mean(log(YGMPE)));
    
for i=1:length(GMPET);
if GMPET(i) <= 0.02
    GMPEStdln(i) = 0.24;
elseif GMPET(i) <=2
    GMPEStdln(i) = 0.08*log10(GMPET(i))+0.376;
else
    GMPEStdln(i) = 0.4;
end
end
    
UB = exp(log(GMPEAve)+GMPEStdln);
LB = exp(log(GMPEAve)-GMPEStdln);    
    
colorAVE = [0 0 0];
FS = 14;

%% Number of rows
%%[nr, np] = size(AS08);
[nr, np] = size(GMPEdata(:,:,1));

%% Plot All GMPE relation, mean, range, and boxplot for simulation results.
    
figure1=figure;
    
ll(1) = plot(GMPET,YGMPE(1,:),'linewidth',1.8,'color',RGB_xls(1,:)); hold all
for i=2:NumGMPE;
ll(i) = plot(GMPET,YGMPE(i,:),'linewidth',1.8,'color',RGB_xls(i,:)); hold all
end
ll(i+1) = plot(GMPET, GMPEAve,'linewidth',2,'color',colorAVE);

h=boxplot(Ysim,'positions',boxposition,'widths',boxwidth,'whisker',6);
hm=findobj('Tag','Median');
hb=findobj('Tag','Box');
ho=findobj('Tag','Outliers');
hw=findobj('Tag','Upper Whisker','-or','Tag','Lower Whisker');
set(hb,'LineStyle','-','LineWidth',1.5);
set(hw,'LineStyle','-','LineWidth',1.1);
set(hm,'Visible','off');
set(ho,'Visible','off');
scatter(X,meansim,50,'r','s','filled'); %plot medians again to make larger
    
GMPELabels = GMPELabels';
GMPELabels(end+1) = {'Mean of GMPE Models'};

legend(ll,GMPELabels,...
       'FontSize',FS,'location','SouthWest'); 

set(gca,'YScale','log','YMinorTick','on',...
    'YMinorGrid','on',...
    'XTickLabel',{'0.01','0.1','1','10'},...
    'XTick',[0.01 0.1 1 10],...
    'XScale','log',...
    'XMinorTick','on',...
    'XMinorGrid','on',...
    'FontSize',14);
xlim([0.01,20]);ylim([ymin ymax]); box on; grid on;

% Create xlabel
xlabel('Period (sec)','Units','points','FontSize',FS);

% Create ylabel
ylabel('PSA (g)','FontSize',FS);

% Create title
title(PlotTitle,'FontWeight','bold','FontSize',FS)

Plot_title = [ OUTFile '_w_gmpe.pdf'];
print(figure1,'-dpdf', '-painters', '-r600', Plot_title)
close

%% Plot GMPE points and range overlain by boxplot for simulations.
    
figure2=figure;
for i=1:nr;
  dd(1) = plot(GMPET,GMPEdata(i,:,1),PLOT_markers(1),'linewidth',2,'markerfacecolor',RGB_xls(1,:),'color',RGB_xls(1,:)); hold all
end
for i=1:nr;
for j=2:NumGMPE;
  dd(j) = plot(GMPET,GMPEdata(i,:,j),PLOT_markers(j),'linewidth',2,'markerfacecolor',RGB_xls(j,:),'color',RGB_xls(j,:));
end
end
dd(j+1) = plot(GMPET,GMPEAve,'linewidth',2,'color',colorAVE);

h=boxplot(Ysim,'positions',boxposition,'widths',boxwidth,'whisker',6);
hm=findobj('Tag','Median');
hb=findobj('Tag','Box');
ho=findobj('Tag','Outliers');
hw=findobj('Tag','Upper Whisker','-or','Tag','Lower Whisker');
set(hb,'LineStyle','-','LineWidth',1.5);
set(hw,'LineStyle','-','LineWidth',1.1);
set(hm,'Visible','off');
set(ho,'Visible','off');
scatter(X,meansim,50,'r','s','filled'); %plot medians again to make larger
    
legend(dd,GMPELabels,...
       'FontSize',FS,'location','SouthWest'); 

set(gca,'YScale','log','YMinorTick','on',...
    'YMinorGrid','on',...
    'XTickLabel',{'0.01','0.1','1','10'},...
    'XTick',[0.01 0.1 1 10],...
    'XScale','log',...
    'XMinorTick','on',...
    'XMinorGrid','on',...
    'FontSize',14);

xlim([0.01,20]);ylim([ymin ymax]); box on; grid on;

% Create xlabel
xlabel('Period (sec)','Units','points','FontSize',FS);

% Create ylabel
ylabel('PSA (g)','FontSize',FS);

% Create title
title(PlotTitle,'FontWeight','bold','FontSize',FS)

Plot_title = [ OUTFile '_w_dots.pdf'];
print(figure2,'-dpdf', '-painters', '-r600', Plot_title)
close

exit
