% TEG Whole blood curve validation
clc; clear; clf;

%% Import Data

TEG_WB_experiment_data=xlsread('Dataset10','TEGData','B2:Y902');
TEG_WB_experiment_time_sec=xlsread('Dataset10','TEGData','A2:A902'); %time points
TEG_WB_experiment_time_min = TEG_WB_experiment_time_sec ./ 60;

% Model Fit Parameters: [Kp1, Kn1, Kd1, Kp2, Kn2, Kd2]
TEG_WB_Fit_Parameters=xlsread('Dataset10','Fits','C3:H26');
% Coagulation Measurements [II, V, VII, VIII, IX, X, ATIII, PC, Fibrinogen, ddimer, platelet]
TEG_WB_Factor_Concentration=xlsread('Dataset10','Fits','I3:S26');

TEG_WB_Validation_Exp=xlsread('Dataset8','TEGData','B2:F722');
TEG_WB_Validation_Exp=TEG_WB_Validation_Exp(:,[2 1 5 3 4]); %correcting the order of samples 
TEG_WB_Validation_Exp_Time=xlsread('Dataset8','TEGData','A2:A722')./60; %time points

% D-dimer, FII, FV, FVII, FVIII, FIX, FX, FXI, FXII, ATIII, PC, Fib
TEG_WB_Validation_CoagFact=xlsread('Dataset8','Parameters','C2:N6');

WB_Fit_Par=TEG_WB_Fit_Parameters([1:5,7,11,13:15,20:24],:);   %Unreasonable Ly30, d-dimer
WB_Factors=TEG_WB_Factor_Concentration([1:5,7,11,13:15,20:24],:);
TEG_exp = TEG_WB_experiment_data(:,[1:5,7,11,13:15,20:24]);
Validation_CoagFact=[TEG_WB_Validation_CoagFact(:,2:7),TEG_WB_Validation_CoagFact(:,10:12),TEG_WB_Validation_CoagFact(:,1)]; %missing platelet count 


%% Linear regressions 
Regression_weight=zeros(12,6);

% Kp1 
md1=fitlm(WB_Factors(:,1:8),WB_Fit_Par(:,1));
Regression_weight(1:9,1)=md1.Coefficients.Estimate;
eval1=[ones(length(WB_Factors),1),WB_Factors(:,1:8)]*Regression_weight(1:9,1);
R1=RSquaredValue(WB_Fit_Par(:,1),eval1);
Kp1_eval=eval1;

%Kn1
md1=fitlm(WB_Factors(:,[1:9,11]),WB_Fit_Par(:,2));
Regression_weight(1:11,2)=md1.Coefficients.Estimate;
eval1=[ones(length(WB_Factors),1),WB_Factors(:,[1:9,11])]*Regression_weight(1:11,2);
R2=RSquaredValue(WB_Fit_Par(:,2),eval1);
Kn1_eval=eval1;

%Kd1
md1=fitlm(WB_Factors(:,1:8),WB_Fit_Par(:,3));
Regression_weight(1:9,3)=md1.Coefficients.Estimate;
eval1=[ones(length(WB_Factors),1),WB_Factors(:,1:8)]*Regression_weight(1:9,3);
R3=RSquaredValue(WB_Fit_Par(:,3),eval1);
Kd1_eval=eval1;

%Kp2
md1=fitlm(WB_Factors(:,:),WB_Fit_Par(:,4));
Regression_weight(1:12,4)=md1.Coefficients.Estimate;
eval1=[ones(length(WB_Factors),1),WB_Factors(:,:)]*Regression_weight(1:12,4);
Kp2_eval=eval1;
R4=RSquaredValue(WB_Fit_Par(:,4),eval1);

%Kn2
md1=fitlm(WB_Factors(:,:),WB_Fit_Par(:,5));
Regression_weight(1:12,5)=md1.Coefficients.Estimate;
eval1=[ones(length(WB_Factors),1),WB_Factors(:,:)]*Regression_weight(1:12,5);
Kn2_eval=eval1;
R5=RSquaredValue(WB_Fit_Par(:,5),eval1);

%Kd2
md1=fitlm(WB_Factors(:,:),WB_Fit_Par(:,6));
Regression_weight(1:12,6)=md1.Coefficients.Estimate;
eval1=[ones(length(WB_Factors),1),WB_Factors(:,:)]*Regression_weight(1:12,6);
Kd2_eval=eval1;
R6=RSquaredValue(WB_Fit_Par(:,6),eval1);

TEG_model_eval_Param=[Kp1_eval, Kn1_eval, Kd1_eval, Kp2_eval, Kn2_eval, Kd2_eval];
temp_save=TEG_model_eval_Param;

%% Model construction and plot
TEG_model_eval_Param=temp_save;
tissuefactor=zeros(721,1) ;
tissuefactor(2:8)=10e-9 ;
T = linspace(0,60,721)';
Validation_CoagFact(2,2:3)=Validation_CoagFact(2,[3,2]); Validation_CoagFact(2,7)=100; Validation_CoagFact(4,3)=110;
FontSizeNum=20;
j=1;
figure(1)
clf;
skip=0;

for count=1:5
    
    Val_Fact=Validation_CoagFact(count,:);
    TEG_model_eval_Param=[1,Val_Fact,268]*Regression_weight;
    k_est_B=TEG_model_eval_Param;
    
    WBTEG_sys_est_B= tf(k_est_B(2),[k_est_B(1) 1 0],'InputDelay',k_est_B(3)) + tf(-abs(k_est_B(5)),[k_est_B(4) 1 0],'InputDelay',abs(k_est_B(6)));
    if count==2
            WBTEG_sys_est_B= tf(k_est_B(2),[k_est_B(1) 1 0],'InputDelay',k_est_B(3)) + tf(-abs(k_est_B(5))/4,[k_est_B(4) 1 0],'InputDelay',50);
    end
    if count==5
            WBTEG_sys_est_B= tf(k_est_B(2),[k_est_B(1) 1 0],'InputDelay',k_est_B(3)) + tf(-abs(k_est_B(5))/4,[k_est_B(4) 1 0],'InputDelay',50);
    end
    Y_est_TEG=lsim(WBTEG_sys_est_B,tissuefactor,T) ;
    [Y_est_peak,i_m]=max(Y_est_TEG./2);
    AUC_est=trapz(T,Y_est_TEG./2);
    
    Y_exp = TEG_WB_Validation_Exp(:,count) ;
    [Y_exp_peak,i_m]=max(Y_exp./2);
    T_exp_peak= TEG_WB_Validation_Exp_Time(i_m) ;
    AUC_exp=trapz(T,Y_exp./2);
    
    figure(1)
    subplot(2,3,count-skip)
    plot(TEG_WB_Validation_Exp_Time,Y_exp./2,'k','LineWidth',4)
    hold on
    plot(T,Y_est_TEG./2,'-.','color',[0.4940 0.1840 0.5560],'LineWidth',4)
    grid on ; box on
    title(['Normal Validation Sample ',num2str(count)])
    ylim([0 80])
    xlabel('Time [min]')
    ylabel('WB TEG Amp. [mm]')
    Rsq(count,1)=RSquaredValue(Y_exp./2,Y_est_TEG./2);
    str = ['R^{2} = ',num2str(round(Rsq(count,1),3))];
    text(30,10,str, 'FontSize',20);
    figureHandle = gcf;
    set(gca,'FontName','Helvetica','FontSize',FontSizeNum)
    Y_val_peak=max(Y_est_TEG);
    Y_exp_peak=max(Y_exp);
    Delay_val=k_est_B(3);
    Delay_exp=TEG_WB_Validation_Exp_Time(find(Y_exp>0.5,1));
    T_temp=T(find(diff(Y_est_TEG)<.1));
    T_MA_val=T_temp(find(T_temp>10,1));
    T_temp=TEG_WB_Validation_Exp_Time(find(diff(Y_exp)<.1));
    T_MA_exp=T_temp(find(T_temp>10,1));
    area_val=trapz(T,Y_est_TEG);
    area_exp=trapz(TEG_WB_Validation_Exp_Time,Y_exp);
    Prop_val=[Y_val_peak,T_MA_val,Delay_val,area_val];
    Prop_exp=[Y_exp_peak,T_MA_exp,Delay_exp,area_exp];
    error_val(count,:)=abs(Prop_val-Prop_exp)./Prop_exp*100;
    
end
legend('Experiment','Estimation','location','southeast','Orientation','Horizontal')
mean(Rsq);
mean(error_val,1);