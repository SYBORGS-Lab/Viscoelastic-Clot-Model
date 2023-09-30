%Rapid TEG and CN TEG Comparison Group
clc; clear; clf;
%% Import Data

% Citrated Rapid: R Time, Angle, MA, LY30
% Citrated Functional Fibrinogen: R Time, Angle, MA, TMA, LY30, CL30, FLEV
% Citrated Native: R Time, K Time, Angle, MA, TMA, LY30, CL30
[COMBAT_TEG_Exp, Header_TEG_Exp]=xlsread('Dataset2','Dataset6CN','B1:S102'); 
%Platelet count
[LAB_Platelet, ~]=xlsread('Dataset2','Dataset6CN','U4:U102'); 

% Divide up the data
DATA_FF_CN=COMBAT_TEG_Exp(:,[5:18]);
DATA_FF=DATA_FF_CN(:,1:7);
DATA_CN=DATA_FF_CN(:,8:14);
DATA_RapidTEG=COMBAT_TEG_Exp(:,1:4);

%% Rapid TEG
FigLineSize=3;
FontSizeNum=30;
MarkSIze=10;

%Fig 1 - Max Amplitude
figure(3)
clf
xDataVal=DATA_CN(:,4);
yDataVal=DATA_RapidTEG(:,3);  
subplot(2,2,1)
plot(xDataVal,yDataVal,'k^','MarkerSize',MarkSIze,'MarkerFaceColor',[0 0.4470 0.7410])
xlabel('CN TEG'); ylabel('Rapid TEG'); title('Maximum Amplitude')

% %linear fit with no y intercept
lin_func_eq=@(a,x) a(1).*x;
a0=1;

%y=x fit
ylin = xDataVal;
hold on;
plot(xDataVal,ylin,'r')
R2_MA=RSquaredValue(yDataVal,ylin)
R1_MA=sqrt(R2_MA);
str=['y   = x'];
text(22.4,75,str,'FontSize',20);
str=['R^{2} = ', num2str(round(R2_MA,3))];
text(22,68,str,'FontSize',20);
grid on; box on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)


%fig 2 - R Time
xDataVal=DATA_CN(:,1);
yDataVal=DATA_CN(:,1)./DATA_RapidTEG(:,1); 
subplot(2,2,2)
RD=find(xDataVal>14);
yDataVal(RD)=[];
xDataVal(RD)=[];
figure(3)
plot(xDataVal,yDataVal,'kd','MarkerSize',MarkSIze,'MarkerFaceColor',[0.8500 0.3250 0.0980])
xlabel('CN TEG'); ylabel('CN TEG/Rapid TEG'); title('R Time')
c1=polyfit(xDataVal,yDataVal,1);
y1=polyval(c1,xDataVal);
hold on; plot(xDataVal,y1,'r')
R2_R=RSquaredValue(yDataVal,y1)
R1_R=sqrt(R2_R);
str=['y   = ',num2str(round(c1(1),3)),'x + ' num2str(round(c1(2),3))];
text(0.57,23,str,'FontSize',20);
str=['R^{2} = ', num2str(round(R2_R,3))];
text(0.5,20,str,'FontSize',20);
ylim([0 25])
xlim([0 14])
grid on; box on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

%fig 3 - LY30
xDataVal=DATA_CN(:,6);
yDataVal=DATA_RapidTEG(:,4); 
subplot(2,2,3)
plot(xDataVal,yDataVal,'ko','MarkerSize',MarkSIze,'MarkerFaceColor',[0.4660 0.6740 0.1880])
xlabel('CN TEG'); ylabel('Rapid TEG'); title('Ly30') 
y1=xDataVal;
hold on; plot(xDataVal,y1,'r')
R2_Ly30=RSquaredValue(yDataVal,y1)
R1_Ly30=sqrt(R2_Ly30);
str=['y   = x'];
text(4.5,92,str,'FontSize',20);
str=['R^{2} = ', num2str(round(R2_Ly30,3))];
text(4,80,str,'FontSize',20);
ylim([0 100])
xlim([0 100])
grid on; box on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

%fig 4 - Alpha Angle
subplot(2,2,4)
xDataVal=DATA_CN(:,3);
yDataVal=DATA_RapidTEG(:,2);  
RD=find(yDataVal>60 & xDataVal<40);
yDataVal(RD)=[];
xDataVal(RD)=[];
plot(xDataVal,yDataVal,'ks','MarkerSize',MarkSIze,'MarkerFaceColor',[0.4940 0.1840 0.5560])
xlabel('CN TEG'); ylabel('Rapid TEG'); title('\alpha Angle') 

%linear fit
lin_func_eq=@(a,x) a(1).*x+a(2);
a0=[1,1];
lin_func_val=lsqcurvefit(lin_func_eq,a0,xDataVal,yDataVal);
xlin=sort(xDataVal);
ylin=lin_func_eq(lin_func_val,xlin);
hold on;
plot(xlin,ylin,'r')
R2_angle=RSquaredValue(yDataVal,lin_func_eq(lin_func_val,xDataVal))
R1_angle=sqrt(R2_angle);
str=['y   = ',num2str(round(lin_func_val(1),3)),'x + ',num2str(round(lin_func_val(2),3))];
text(22.4,79,str,'FontSize',20);
str=['R^{2} = ', num2str(round(R2_angle,3))];
text(22,71,str,'FontSize',20);
grid on; box on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)