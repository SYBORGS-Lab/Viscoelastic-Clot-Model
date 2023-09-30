%% Functional Fibrinogen Plots
clc; clear; clf;

%% Import Data

% Citrated Functional Fibrinogen: R Time, Angle, MA, TMA, LY30, CL30, FLEV
[COMBAT_TEG_Exp, Header_TEG_Exp]=xlsread('Dataset2','Dataset6FF','F1:L66'); 
% Fibrinogen Level
FIBRINOGEN_LEVEL=xlsread('Dataset2','Dataset6FF','T4:T66'); 

DATA_FF_FIB=[COMBAT_TEG_Exp, FIBRINOGEN_LEVEL];
%% Plots
MarkSIze=10;
FontSizeNum=30;

%FLEV vs G_f
figure(1)
subplot(1,2,1)
xAxisData=DATA_FF_FIB(:,3)*5000./(100-DATA_FF_FIB(:,3));
yAxisData=DATA_FF_FIB(:,7);
plot(xAxisData,yAxisData,'k^','MarkerSize',MarkSIze,'MarkerFaceColor',[0 0.4470 0.7410])
c1=polyfit(xAxisData,yAxisData,1);
y1=polyval(c1,xAxisData);
hold on; plot(xAxisData,y1,'r')
R2_FLEV=RSquaredValue(yAxisData,y1);
xlabel('G_f'); ylabel('TEG FLEV');
ylim([0 800]);
str=['y   = ',num2str(round(c1(1),3)),'x + ',num2str(round(c1(2),3))];
text(220,750,str,'FontSize',20);
str=['R^{2} = 0.990'];
text(200,720,str,'FontSize',20);
box on; grid on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

%Fibrinogen vs G_f
figure(1)
subplot(1,2,2)
xAxisData=DATA_FF_FIB(:,3)*5000./(100-DATA_FF_FIB(:,3));
yAxisData=DATA_FF_FIB(:,8);
xDataVal=xAxisData;
yDataVal=yAxisData; 
hold on
plot(xDataVal,yDataVal,'kd','MarkerSize',MarkSIze,'MarkerFaceColor',[0.8500 0.3250 0.0980])
xlabel('G_f'); ylabel('Fibrinogen Level');
ylim([0 400]);
c1=polyfit(xDataVal,yDataVal,1);
y1=polyval(c1,xDataVal);
hold on; plot(xDataVal,y1,'r')
R2_AvFF=RSquaredValue(yDataVal,y1);
R1_AvFF=sqrt(R2_AvFF);
str=['y   = 0.120x + ',num2str(round(c1(2),3))];
text(220,375,str,'FontSize',20);
str=['R^{2} = ',num2str(round(R2_AvFF,3))];
text(200,360,str,'FontSize',20);
box on; grid on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)