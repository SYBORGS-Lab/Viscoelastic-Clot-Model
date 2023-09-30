%% Plot of Fits versus TEG Parameters
clc; clear; clf;

%% Import Data 

% Model Fit Parameters: [Kp1, Kn1, Kd1, Kp2, Kn2, Kd2]
TEG_WB_6_parameters=xlsread('Dataset10','Fits','C3:H26');

% Coagulation Measurements [II, V, VII, VIII, IX, X, ATIII, PC, Fibrinogen, ddimer, platelet]
TEG_24tra_CoagFac_data=xlsread('Dataset10','Fits','I3:S26');

% Machine read TEG Parameters: [R, K, alpha, MA, Ly30]
TEG_24tra_MachProp=xlsread('Dataset10','Fits','U3:Z26');

%% Fits vs Prop
FigLineSize=3;
FontSizeNum=30;
MarkSIze=10;

R1 = zeros(4,1);
R2 = zeros(4,1);

% Kn1 vs MA
figure(1); 
subplot(2,2,1)
xAxisData=TEG_WB_6_parameters(:,2);
yAxisData=TEG_24tra_MachProp(:,4);
plot(xAxisData,yAxisData,'k^','MarkerSize',MarkSIze,'MarkerFaceColor',[0 0.4470 0.7410])
xlabel('K_{n1}'); ylabel('MA');
hold on;
lin_eq=@(b,x) b(1).*x+b(2);
b0=[0,1];
lin_val=lsqcurvefit(lin_eq,b0,xAxisData,yAxisData);
xlin=sort(xAxisData);
ylin=lin_eq(lin_val,xlin);
R2(1)=RSquaredValue(yAxisData,lin_eq(lin_val,xAxisData));
R1(1)=sqrt(R2(1));
plot(xlin,ylin,'r')
str=['y   = 3.006x10^{-9}x - ',num2str(abs(round(lin_val(2),3)))];
text(18450000000,77,str,'FontSize',20);
str=['R^{2} = ',num2str(round(R2(1),3))];
text(18400000000,73.5,str,'FontSize',20);
grid on; box on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)


% Kd1 vs R time
figure(1); 
subplot(2,2,2)
xAxisData=TEG_WB_6_parameters(:,3);
yAxisData=TEG_24tra_MachProp(:,1);
plot(xAxisData,yAxisData,'kd','MarkerSize',MarkSIze,'MarkerFaceColor',[0.8500 0.3250 0.0980])
xlabel('K_{d1}'); ylabel('R Time');
c1=polyfit(xAxisData,yAxisData,1);
y1=polyval(c1,xAxisData);
hold on; plot(xAxisData,y1,'r')
R2(2)=RSquaredValue(yAxisData,y1);
R1(2)=sqrt(R2(2));
str=['y   = ',num2str(round(c1(1),3)),'x + ',num2str(round(c1(2),3))];
text(2.225,5.9,str,'FontSize',20);
str=['R^{2} = ',num2str(round(R2(2),3))];
text(2.2,5.5,str,'FontSize',20);
grid on; box on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)


% Kp1 vs alpha
figure(1); 
subplot(2,2,3) 
xAxisData=TEG_WB_6_parameters(:,1);
yAxisData=TEG_24tra_MachProp(:,3);
RD_1=find(yAxisData<60);
RD=RD_1;
xAxisData(RD)=[];
yAxisData(RD)=[];
plot(xAxisData,yAxisData,'ko','MarkerSize',MarkSIze,'MarkerFaceColor',[0.4660 0.6740 0.1880])
xlabel('K_{p1}'); ylabel('\alpha Angle');
c1=polyfit(xAxisData,yAxisData,1);
y1=polyval(c1,xAxisData);
hold on; plot(xAxisData,y1,'r')
R2(3)=RSquaredValue(yAxisData,y1);
R1(3)=sqrt(R2(3));
str=['y = ',num2str(round(c1(1),3)),'x + ',num2str(round(c1(2),3))];
text(5.1,78,str,'FontSize',20);
str=['R^{2} = ',num2str(round(R2(3),3))];
text(6.3,75.5,str,'FontSize',20);
grid on; box on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

% Degredation (Area under the curve of the lysis component)
tissuefactor=zeros(901,1) ;
tissuefactor(2:8)=10e-9 ;
T = linspace(0,75,901)';
degredation=zeros(24,1);
for cnt=1:24
    sys_test=tf(TEG_WB_6_parameters(cnt,5),[TEG_WB_6_parameters(cnt,4),1,0],'InputDelay',TEG_WB_6_parameters(cnt,6));
    Y_est = lsim(sys_test,tissuefactor,T) ;
    degredation(cnt,1)=trapz(T,Y_est);
end

% Degredation vs LY30
figure(1); 
subplot(2,2,4)
xAxisData=degredation;
yAxisData=TEG_24tra_MachProp(:,5);
plot(xAxisData,yAxisData,'ks','MarkerSize',MarkSIze,'MarkerFaceColor',[0.4940 0.1840 0.5560])
ylabel('Ly30'); xlabel('AUC');
c1=polyfit(xAxisData,yAxisData,1);
y1=polyval(c1,xAxisData);
hold on; plot(xAxisData,y1,'r')
R2(4)=RSquaredValue(yAxisData,y1);
R1(4)=sqrt(R2(4));
str=['y   = 6.000x10^{-3}x - ',num2str(abs(round(c1(2),3)))];
text(110,11,str,'FontSize',20);
str=['R^{2} = ',num2str(round(R2(4),3))];
text(100,9.3,str,'FontSize',20);
grid on; box on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

% Output Table for the R squared values
Plots = ["Kn1 vs MA";"Kd1 vs R time";"Kp1 vs alpha";"Degredation vs LY30"];
rvalues = table(Plots,R2,R1)