% TEG Plasma curve estimation 
clc; clear; clf;

%% Import Data
[TEG_Fits_Pole_SD] = xlsread('Dataset7','Fits','A1:D10');
[SD_Factor] = xlsread('Dataset7','Factor','A1:G10');
[CATData_SD] = xlsread('Dataset7','CATData','B2:K180');
TimeM = xlsread('Dataset7','CATData','A2:A180');
[TEGData_SD] = xlsread('Dataset7','TEGData','B2:K180');

%% Model construction and plot : Fits , Model train on all SD and apply on all SD 

for k=1:3 
    
    %NON-GREEDY Method
    param=TEG_Fits_Pole_SD(:,k);
    fit_param_sorted=pinv(SD_Factor)*param ;
    const_NGree(1,k)=0 ;
    factor_coeff_NGree(:,k)=fit_param_sorted ;
    
        
    % GREEDY Method
    param=TEG_Fits_Pole_SD(:,k);
    [res,fit_param_factor_chosen,i_chosen, fit_param_sorted] = GreedyLSQ(SD_Factor,param) ;
    const_G(1,k)=fit_param_factor_chosen(1) ;
    factor_coeff_G(:,k)=fit_param_sorted ;
    
    % Greedy restrict 
    param= TEG_Fits_Pole_SD(:,k);
    if k==1
        [res,fit_param_factor_chosen,i_chosen, fit_param_sorted] = GreedyLSQ(SD_Factor(:,7),param) ;
    elseif k==2
        [res,fit_param_factor_chosen,i_chosen, fit_param_sorted] = GreedyLSQ(SD_Factor(:,1),param) ;
    else 
        [res,fit_param_factor_chosen,i_chosen, fit_param_sorted] = GreedyLSQ(SD_Factor,param) ;
    end
    const_RsG(1,k)=fit_param_factor_chosen(1) ;
    factor_coeff_RsG(1:7,k)=fit_param_sorted ;
    
end
%% estimation 10 normal 

FontSizeNum=15;

for i=1:10
   
    CAT_exp=CATData_SD(:,i) ;
    
    factor_concen=SD_Factor(i,:);  
    
    ks_fit=TEG_Fits_Pole_SD(i,:) ;
    ks_est_NG=factor_concen*factor_coeff_NGree+const_NGree;
    ks_est_G=factor_concen*factor_coeff_G+const_G;

    ks_est_RsG(1)=SD_Factor(i,1)*factor_coeff_RsG(1,1)+const_RsG(1);
    ks_est_RsG(2)=SD_Factor(i,7)*factor_coeff_RsG(1,2)+const_RsG(2);
    ks_est_RsG(3)=SD_Factor(i,:)*factor_coeff_RsG(:,3)+const_RsG(3);
    
    sys_fit= tf(ks_fit(1),[ks_fit(3) 1 0],'InputDelay',ks_fit(2));
    sys_est_NG= tf(ks_est_NG(1),[ks_est_NG(3) 1 0],'InputDelay',ks_est_NG(2));
    sys_est_G= tf(ks_est_G(1),[ks_est_G(3) 1 0],'InputDelay',ks_est_G(2));
    sys_est_RsG= tf(ks_est_RsG(1),[ks_est_G(3) 1 0],'InputDelay',ks_est_G(2));
    
    T = linspace(0,60,179)';
       
    Y_exp = TEGData_SD(:,i) ;
    [Y_exp_peak,i_m]=max(Y_exp);
    T_exp_peak= TimeM(i_m) ;
    
    Y_fit=lsim(sys_fit,CAT_exp,T) ;
    [Y_fit_peak,i_m]=max(Y_fit);
    T_fit_peak= T(i_m) ;
    
    Y_est_NG=lsim(sys_est_NG,CAT_exp,T) ;
    [Y_est_NG_peak,i_m]=max(Y_est_NG);
    T_est_NG_peak= T(i_m) ;
    
    Y_est_G=lsim(sys_est_G,CAT_exp,T) ;
    [Y_est_G_peak,i_m]=max(Y_est_G);
    T_est_G_peak= T(i_m) ;
    
    Y_est_RsG=lsim(sys_est_RsG,CAT_exp,T) ;
    [Y_est_RsG_peak,i_m]=max(Y_est_RsG);
    T_est_RsG_peak= T(i_m) ;
    
    Y_peak(i,1)= Y_exp_peak ;
    Y_peak(i,2)= Y_fit_peak ;
    Y_peak(i,3)= Y_est_NG_peak ;
    Y_peak(i,4)= Y_est_G_peak ;
    Y_peak(i,5)= Y_est_RsG_peak ;
    
    T_peak(i,1)= T_exp_peak ;
    T_peak(i,2)= T_fit_peak ;
    T_peak(i,3)= T_est_NG_peak ;
    T_peak(i,4)= T_est_G_peak ;
    T_peak(i,5)= T_est_RsG_peak ;
    
    T_delay(i,1)= DetermineDelay(TimeM, Y_exp) ;
    T_delay(i,2)= DetermineDelay(T, Y_fit) ;
    T_delay(i,3)= DetermineDelay(T, Y_est_NG) ;
    T_delay(i,4)= DetermineDelay(T, Y_est_G) ;
    T_delay(i,5)= DetermineDelay(T, Y_est_RsG) ;
    
    R(i,1)=RSquaredValue(Y_exp, Y_fit) ;
    R(i,2)=RSquaredValue(Y_exp, Y_est_NG) ;
    R(i,3)=RSquaredValue(Y_exp, Y_est_G) ;
    R(i,4)=RSquaredValue(Y_exp, Y_est_RsG) ;
 
%Exp vs Estimate 
    figure (2) 
    subplot(2,5,i)
    plot(TimeM,Y_exp,'k',T,Y_est_G,'LineWidth',4)
    grid on; box on
    ylim([0 80])
    str = ['R^{2} = ',num2str(round(R(i,3),3))];
    text(5,75,str);
    str = ['Plasma Sample ',num2str(i)];
    title(str)
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)
    set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

end
legend('Experiment','Estimation','Location','northwest')
subplot(2,5,8); xlabel('Time [min]')
subplot(2,5,1);  ylabel('Amplitude [mm]')
subplot(2,5,6);  ylabel('Amplitude [mm]')
set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)
figure (2) 
set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

%% Validation using 5 normal 

TEG_Plasma_Validation_Exp=xlsread('Dataset8','TEGData','H2:L722');
TEG_Plasma_Validation_Exp=TEG_Plasma_Validation_Exp(:,[2 1 5 3 4]); %correcting the order of samples 
TEG_Plasma_Validation_Exp_Time=xlsread('Dataset8','TEGData','A2:A722')./60; %time points
% D-dimer, FII, FV, FVII, FVIII, FIX, FX, FXI, FXII, ATIII, PC, Fib
TEG_Plasma_Validation_CoagFact=xlsread('Dataset8','Parameters','C2:N6');

Validation_CoagFact=[TEG_Plasma_Validation_CoagFact(:,2:7),TEG_Plasma_Validation_CoagFact(:,10)]; 

% Linear regressions 
Regression_weight=zeros(8,3);

% Kp1 
md1=fitlm(SD_Factor(:,1:7),TEG_Fits_Pole_SD(:,3));
Regression_weight(:,3)=md1.Coefficients.Estimate;
eval1=[ones(length(SD_Factor),1),SD_Factor(:,1:7)]*Regression_weight(:,3);
R1=RSquaredValue(TEG_Fits_Pole_SD(:,3),eval1);
Kp1_eval=eval1;

%Kn1
md1=fitlm(SD_Factor(:,1:7),TEG_Fits_Pole_SD(:,1));
Regression_weight(:,1)=md1.Coefficients.Estimate;
eval1=[ones(length(SD_Factor),1),SD_Factor(:,1:7)]*Regression_weight(:,1);
R2=RSquaredValue(TEG_Fits_Pole_SD(:,1),eval1);
Kn1_eval=eval1;

%Kd1
md1=fitlm(SD_Factor(:,1:7),TEG_Fits_Pole_SD(:,2));
Regression_weight(:,2)=md1.Coefficients.Estimate;
eval1=[ones(length(SD_Factor),1),SD_Factor(:,1:7)]*Regression_weight(:,2);
R3=RSquaredValue(TEG_Fits_Pole_SD(:,2),eval1);
Kd1_eval=eval1;


TEG_Plasmamodel_eval_Param=[Kn1_eval, Kd1_eval, Kp1_eval];
temp_save=TEG_Plasmamodel_eval_Param;


CAT_Normal_5validation=xlsread('Dataset8','TEGData','N2:R125')./1;
CAT_Normal_5validation=[CAT_Normal_5validation;zeros(57,5)];
Validation_CoagFact=[TEG_Plasma_Validation_CoagFact(:,2:7),TEG_Plasma_Validation_CoagFact(:,10)]; 
Validation_CoagFact(2,2:3)=Validation_CoagFact(2,[3,2]); Validation_CoagFact(2,7)=100; Validation_CoagFact(4,3)=110;
FontSizeNum=20;
figure(8)
clf
clear error_val

for count=1:5
    CAT_exp=CAT_Normal_5validation(:,count);
    Val_Fact=Validation_CoagFact(count,:);
    K_plasma_Val=[1,Val_Fact]*Regression_weight;
    
    sys_val_est= tf(K_plasma_Val(1)/1.2,[abs(K_plasma_Val(3))/12 1 0],'InputDelay',min(10,K_plasma_Val(2))-3);
    T = linspace(0,60,181)';
    Y_val_est=lsim(sys_val_est,CAT_exp,T);
    
    Y_val_peak=Y_val_est(91);
    Y_exp_peak=TEG_Plasma_Validation_Exp(361,count);
    Delay_val=min(10,K_plasma_Val(2))-3;
    Delay_exp=TEG_Plasma_Validation_Exp_Time(find(TEG_Plasma_Validation_Exp(:,count)>0.5,1));
    T_temp=T(find(diff(Y_val_est)<.1));
    T_MA_val=T_temp(find(T_temp>10,1));
    T_temp=TEG_Plasma_Validation_Exp_Time(find(diff(TEG_Plasma_Validation_Exp(:,count))<.1));
    T_MA_exp=T_temp(find(T_temp>10,1));
    area_val=trapz(T,Y_val_est);
    area_exp=trapz(TEG_Plasma_Validation_Exp_Time,TEG_Plasma_Validation_Exp(:,count));
    
    Prop_val=[Y_val_peak,T_MA_val,Delay_val,area_val];
    Prop_exp=[Y_exp_peak,T_MA_exp,Delay_exp,area_exp];
    error_val(count,:)=abs(Prop_val-Prop_exp)./Prop_exp*100;
    Prop_Fold_Error_TABLE(count,:)=mean(error_val(count,:),1);
    R_sqr(count)=RSquaredValue(TEG_Plasma_Validation_Exp(1:4:end,count)./2,Y_val_est./2);
    
    figure(8)
    subplot(2,3,count)
    plot(TEG_Plasma_Validation_Exp_Time,TEG_Plasma_Validation_Exp(:,count)./2,'k','LineWidth',3)
    hold on; grid on; box on
    plot(T,Y_val_est./2,'-.','color',[0.4940 0.1840 0.5560],'LineWidth',3)
    ylim([0 50])
    xlabel('Time [min]')
    ylabel('Plasma TEG Amp. [mm]')
    str = ['R^{2} = ',num2str(round(R_sqr(count),3))];
    text(30,5,str, 'FontSize',20);
    figureHandle = gcf;
    set(gca,'FontName','Helvetica','FontSize',FontSizeNum)
    title(sprintf('Normal Validation Sample %d', count))
end
legend('Experiment','Estimation')
MA_5norm_mean=mean(error_val(:,1));
MATime_5norm_mean=mean(error_val(:,2));
Delay_5norm_mean=mean(error_val(:,3));
AUC_5norm_mean=mean(error_val(:,4));
mean(R_sqr);