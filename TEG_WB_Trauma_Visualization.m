% TEG Whole blood fits reconstruction 
clc; clear; clf;
%% Import Data

TEG_WB_experiment_data=xlsread('Dataset10','TEGData','B2:Y902');
TEG_WB_experiment_time_sec=xlsread('Dataset10','TEGData','A2:A902'); %time points
TEG_WB_experiment_time_min = TEG_WB_experiment_time_sec ./ 60;

TEG_exp = TEG_WB_experiment_data(:,[1:2,5:6,8:13,15,18:24]);

%% Sample TEGS 1, 11
figure(1)
plot(TEG_WB_experiment_time_min,TEG_exp(:,1),TEG_WB_experiment_time_min,TEG_exp(:,11),'LineWidth',3)
legend('Trauma sample 1','Trauma sample 2');
ax = gca;
ax.FontSize = 20; 
grid on
box on
xlabel('Time [min]')
ylabel('Amplitude [mm]')