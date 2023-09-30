%% Platelet Mapping
clc; clear; clf;

%% Import Data

% Citrated Functional Fibrinogen: R Time, Angle, MA, TMA, LY30, CL30, FLEV
% Citrated Native: R Time, K Time, Angle, MA, TMA, LY30, CL30
% Platelet Mapping: MA FIB, MA ADP, % ADP Inhibition
[COMBAT_TEG_Exp, Header_TEG_Exp]=xlsread('Dataset2','Dataset11','B1:R59'); 
%Platelet count
[LAB_Platelet, ~]=xlsread('Dataset2','Dataset11','S4:S59'); 

DATA_FF_CN_PLT=[COMBAT_TEG_Exp, LAB_Platelet];

%Data Clean Up - delete where MA FIB is larger than MA ADP
[row0, ~] = find(DATA_FF_CN_PLT(:,16)<=DATA_FF_CN_PLT(:,15));
DATA_FF_CN_PLT(row0,:)=[];


MA_PLT_Thromb=DATA_FF_CN_PLT(:,16).*(1./(1-DATA_FF_CN_PLT(:,17)/100));
%CN MA - FF MA
Platelet_FF_CN=DATA_FF_CN_PLT(:,11)-DATA_FF_CN_PLT(:,3);
%CN MA / FF MA
Platelet_FF_CN_Ratio=DATA_FF_CN_PLT(:,11)./DATA_FF_CN_PLT(:,3);
%MA ADP - MA FIB
Platelet_PLT=DATA_FF_CN_PLT(:,16)-DATA_FF_CN_PLT(:,15);
%MA ADP / MA FIB
Platelet_PLT_Ratio=DATA_FF_CN_PLT(:,16)./DATA_FF_CN_PLT(:,15);
%Just all the lab plt counts
LAB_Platelet_Count=DATA_FF_CN_PLT(:,18);
MA_CN=DATA_FF_CN_PLT(:,11);
MA_FF=DATA_FF_CN_PLT(:,3);
%MA for FIB plt
MA_PLT_A=DATA_FF_CN_PLT(:,15);
LAB_Platelet_Count_inh=LAB_Platelet_Count.*(1-(DATA_FF_CN_PLT(:,17)./100));
INH_Per=DATA_FF_CN_PLT(:,17);
row0=find(Platelet_PLT>mean(Platelet_PLT)+2*std(Platelet_PLT));
row1=find(MA_PLT_Thromb>mean(MA_PLT_Thromb)+2*std(MA_PLT_Thromb));
row0=unique([row0;row1]);
Platelet_FF_CN(row0,:)=[];
Platelet_FF_CN_Ratio(row0,:)=[];
Platelet_PLT(row0,:)=[];
Platelet_PLT_Ratio(row0,:)=[];
LAB_Platelet_Count(row0,:)=[];
LAB_Platelet_Count_inh(row0,:)=[];
MA_PLT_Thromb(row0,:)=[];
MA_CN(row0,:)=[];
MA_FF(row0,:)=[];
MA_PLT_A(row0,:)=[];
INH_Per(row0,:)=[];

% Linear fits set up
lin_eq_no_yint=@(a,x) a(1).*x;
lin_eq=@(b,x) b(1).*x+b(2);
a0=1;
b0=[1,1];

%%  Plots
FontSizeNum=22;
MarkSIze=10;

%Graph A - Platelet vs CN
figure(1); subplot(2,2,1);
xDataVal=MA_CN;
yDataVal=MA_PLT_Thromb;
RD=find(yDataVal<20 | yDataVal>90 | xDataVal<50);
yDataVal(RD)=[];
xDataVal(RD)=[];
plot(xDataVal,yDataVal,'k^','MarkerSize',MarkSIze,'MarkerFaceColor',[0 0.4470 0.7410])
xlabel('Citrated Native'); ylabel('Platelet Mapping'); title('A: Whole Blood MA (Uninhibited)');
hold on;
%y=mx fit
lin_func_val=lsqcurvefit(lin_eq_no_yint,a0,xDataVal,yDataVal);
xlin=sort(xDataVal);
ylin=lin_eq_no_yint(lin_func_val,xlin);
plot(xlin,ylin,'r')
R2_A=RSquaredValue(yDataVal,lin_eq_no_yint(lin_func_val,xDataVal));
R1_A=sqrt(R2_A);
str=['y   = ',num2str(round(lin_func_val(1),3)),'x'];
text(47.2,82,str,'FontSize',20);
str=['R^{2} = ', num2str(round(R2_A,3))];
text(47,78,str,'FontSize',20);
box on; grid on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)
ylim([45 85])
xlim([45 85])


%Graph B - Platelet vs FF
figure(1); subplot(2,2,2);
xAxisData=MA_FF;
yAxisData=MA_PLT_A;
RD=find(yAxisData>15 | yAxisData<0.3378*xAxisData-4);
xAxisData(RD,:)=[];
yAxisData(RD,:)=[];
xDataVal=xAxisData;
yDataVal=yAxisData; 
hold on
plot(xDataVal,yDataVal,'kd','MarkerSize',MarkSIze,'MarkerFaceColor',[0.8500 0.3250 0.0980])
xlabel('Functional Fibrinogen'); ylabel('Platelet Mapping'); title('B: MA from Fibrin');
hold on;
%linear fit - no y intercept
lin_func_val=lsqcurvefit(lin_eq_no_yint,a0,xDataVal,yDataVal);
xlin=sort(xDataVal);
ylin=lin_eq_no_yint(lin_func_val,xlin);
R2_B=RSquaredValue(yDataVal,lin_eq_no_yint(lin_func_val,xDataVal));
R1_B=sqrt(R2_B);
plot(xlin,ylin,'r')
str=['y   = ',num2str(round(lin_func_val(1),3)),'x'];
text(2.2,13.5,str,'FontSize',20);
str=['R^{2} = ', num2str(round(R2_B,3))];
text(2,12,str,'FontSize',20);
box on; grid on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)


%Graph D - Platelet contribution to MA vs Platelet Count
figure(1); subplot(2,2,4);
xAxisData=MA_PLT_Thromb./MA_PLT_A;
yAxisData=(INH_Per./100).*LAB_Platelet_Count;
RD=find(yAxisData<5.0986*xAxisData-30);
xAxisData(RD,:)=[];
yAxisData(RD,:)=[];
xDataVal=yAxisData;
yDataVal=xAxisData; 
hold on;
plot(xDataVal,yDataVal,'ks','MarkerSize',MarkSIze,'MarkerFaceColor',[0.4940 0.1840 0.5560])
xlabel('Uninhibited Platelet Count'); ylabel('Contribution to MA'); title('D: Uninhibited Platelet');
hold on;
%linear fit with no y intercept
lin_func_val=lsqcurvefit(lin_eq_no_yint,a0,xDataVal,yDataVal);
xlin=sort(xDataVal);
ylin=lin_eq_no_yint(lin_func_val,xlin);
R2_D=RSquaredValue(yDataVal,lin_eq_no_yint(lin_func_val,xDataVal));
R1_D=sqrt(R2_D);
plot(xlin,ylin,'r')
str=['y   = ',num2str(round(lin_func_val(1),3)),'x'];
text(16,54.5,str,'FontSize',20);
str=['R^{2} = ', num2str(round(R2_D,3))];
text(15,48.5,str,'FontSize',20);
box on; grid on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)


%Graph C - % MA Reduction vs. % Inhibition
figure(1); subplot(2,2,3);
xAxisData=(INH_Per);
yAxisData=((MA_PLT_Thromb-MA_PLT_A)./(MA_PLT_Thromb))*100;
xDataVal=xAxisData;
yDataVal=yAxisData; 
hold on;
plot(xDataVal,yDataVal,'ko','MarkerSize',MarkSIze,'MarkerFaceColor',[0.4660 0.6740 0.1880])
xlabel('% Inhibition'); ylabel('% MA Reduction'); title('C: ADP Inhibition');
hold on;
%linear fit thru 100,100
lin_eq=@(b,x) (1-(b(1)/100)).*x+b(1);
b0=1;
lin_func_val=lsqcurvefit(lin_eq,b0,xDataVal,yDataVal);
xlin=sort(xDataVal);
ylin=lin_eq(lin_func_val,xlin);
R2_C=RSquaredValue(yDataVal,lin_eq(lin_func_val,xDataVal));
R1_C=sqrt(R2_C);
m= (1-(lin_func_val(1)/100));
disp(lin_eq(lin_func_val,100))
plot(xlin,ylin,'r')
str=['y = ',num2str(round(m,3)),'x + ',num2str(round(lin_func_val(1),3))];
text(56,69,str,'FontSize',20);
str=['R^{2} = ', num2str(round(R2_C,3))];
text(75,65,str,'FontSize',20);
box on; grid on;
figureHandle = gcf;
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)