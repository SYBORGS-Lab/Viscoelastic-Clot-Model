%% 10 Normal Sample Data
clc; clear; clf;

%% Import Data
CATtime = xlsread('Dataset7','CATData','A2:A180');
[CATSample] = xlsread('Dataset7','CATData','B2:K180');
TEGtime = xlsread('Dataset7','TEGData','A2:A180');
[TEGSample] = xlsread('Dataset7','TEGData','B2:K180');
[Fits] = xlsread('Dataset7','Fits','A1:E10');
[Factors] = xlsread('Dataset7','Factor','A1:G10');

%% Plots
FontSizeNum=15;

%Thrombin vs time (CAT)
figure(1)
for i = 1:1:10
    hold on
    plot(CATtime,CATSample(:,i),'LineWidth',4)
end
title('CAT Data - Normal Single Donor Plasma')
xlabel('Time [min]')
ylabel('Thrombin Generation [nM]')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)
hold off

%MA vs time (TEG)
figure(2)
for i = 1:1:10
    hold on
    plot(TEGtime,TEGSample(:,i),'LineWidth',4)
end
title('TEG Data - Normal Single Donor Plasma')
xlabel('Time [min]')
ylabel('TEG Amplitude [mm]')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)
hold off

%knp vs kd
figure(3)
plot(Fits(:,5),Fits(:,1),'linestyle','none','marker','o','MarkerFaceColor','magenta','MarkerEdgeColor','magenta')
hold on
% linear fit
lin_func_eq=@(a,x) a(1).*x+a(2);
a0=[1,1];
% remove outlier
[NewFits] = [Fits(1:2,:); Fits(4:10,:)];
lin_func_val=lsqcurvefit(lin_func_eq,a0,NewFits(:,5),NewFits(:,1));
xlin=sort(NewFits(:,5));
ylin=lin_func_eq(lin_func_val,xlin);
plot(xlin,ylin,'k','LineWidth',2)
hold off
ylim([0 0.05])
xlabel('CAT Model Parameter K_{d}')
ylabel('TEG Model Parameter K_{n_{p}}')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

%knp vs ATIII
figure(4)
plot(Factors(:,7),Fits(:,1),'linestyle','none','marker','o','MarkerFaceColor','blue')
hold on
% linear fit
lin_func_eq=@(a,x) a(1).*x+a(2);
a0=[1,1];
% remove outliers
[NewFac] = [Factors(1:7,:); Factors(10,:)];
[NewFits] = [Fits(1:7,:); Fits(10,:)];
lin_func_val=lsqcurvefit(lin_func_eq,a0,NewFac(:,7),NewFits(:,1));
xlin=sort(NewFac(:,7));
ylin=lin_func_eq(lin_func_val,xlin);
plot(xlin,ylin,'k','LineWidth',2)
hold off
ylim([0 0.05])
xlabel('ATIII')
ylabel('TEG Model Parameter K_{n_{p}}')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)