% TEG Whole blood curve estimation  
clc; clear; clf;

%% Import Data

TEG_WB_experiment_data=xlsread('Dataset10','TEGData','B2:Y902');
TEG_WB_experiment_time_sec=xlsread('Dataset10','TEGData','A2:A902'); %time points
TEG_WB_experiment_time_min = TEG_WB_experiment_time_sec ./ 60;

% Model Fit Parameters: [Kp1, Kn1, Kd1, Kp2, Kn2, Kd2]
TEG_WB_Fit_Parameters=xlsread('Dataset10','Fits','C3:H26');
% Coagulation Measurements [II, V, VII, VIII, IX, X, ATIII, PC, Fibrinogen, ddimer, platelet]
TEG_WB_Factor_Concentration=xlsread('Dataset10','Fits','I3:S26');

WB_Fit_Par=TEG_WB_Fit_Parameters([1:5,7,11,13:15,20:24],:);   %Unreasonable Ly30, d-dimer
WB_Factors=TEG_WB_Factor_Concentration([1:5,7,11,13:15,20:24],:);
TEG_exp = TEG_WB_experiment_data(:,[1:5,7,11,13:15,20:24]);

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
tissuefactor=zeros(901,1) ;
tissuefactor(2:8)=10e-9 ;
T = linspace(0,75,901)';

TEG_model_eval_Param=abs(TEG_model_eval_Param);
TEG_model_eval_Param([5,6],[5])=TEG_model_eval_Param([5,6],[5])*2;
TEG_model_eval_Param([12,13],[4])=TEG_model_eval_Param([12,13],[4])./50;
FontSizeNum=25;
j=1;
figure(1)
clf
skip=0;

for count=1:15 
    k_est_B=TEG_model_eval_Param(count,:);
    WBTEG_sys_est_B= tf(k_est_B(2),[k_est_B(1) 1 0],'InputDelay',k_est_B(3)) + tf(-k_est_B(5),[k_est_B(4) 1 0],'InputDelay',k_est_B(6));
    Y_est_TEG=lsim(WBTEG_sys_est_B,tissuefactor,T) ;
    [Y_est_peak,i_m]=max(Y_est_TEG./2);
    AUC_est=trapz(T,Y_est_TEG./2);
    
    Y_exp = TEG_exp(:,count) ;
    [Y_exp_peak,i_m]=max(Y_exp./2);
    T_exp_peak= TEG_WB_experiment_time_min(i_m) ;
    AUC_exp=trapz(T,Y_exp./2);
    
    subplot(3,5,count-skip)
    plot(TEG_WB_experiment_time_min,Y_exp./2,'k','LineWidth',4)
    hold on
    plot(T,Y_est_TEG./2,'color',[0 0.4470 0.7410],'LineWidth',4)
    grid on ; box on
    title(['Patient ',num2str(count)])
    ylim([0 80])
    xlim([0 75])
    figureHandle = gcf;
    set(gca,'FontName','Helvetica','FontSize',FontSizeNum)
    Rsq(count,1)=RSquaredValue(Y_exp./2,Y_est_TEG./2);
    str=['R^{2} = ', num2str(Rsq(count,1),'%0.3f')];
    text(18,15,str,'FontSize',20);
    
    %Func_out=[TEG_R,TEG_K,TEG_alpha,TEG_MA,TEG_Ly30,TEG_MAtime]
    Func_out=TEG_Graph_Property_Identifier(T,Y_est_TEG,2);
    ParamTable_est(count,:)=[Y_est_peak,AUC_est,Func_out];
    Func_out=TEG_Graph_Property_Identifier(TEG_WB_experiment_time_min,Y_exp,1);
    ParamTable_exp(count,:)=[Y_exp_peak,AUC_exp,Func_out];
end
legend('Experiment','Estimation','location','southeast','Orientation','Horizontal')
subplot(3,5,13); xlabel('Time [min]')
subplot(3,5,6);  ylabel('Amplitude [mm]')
mean(Rsq);

%Error_Table=[MA, AUC, R Time, k, alpha angle, MA, LY30, MA Time]
Error_TABLE=((abs(ParamTable_est-ParamTable_exp))./ParamTable_exp)*100;
mean(Error_TABLE,1);