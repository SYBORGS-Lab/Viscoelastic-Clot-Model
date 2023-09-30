%% Plot of the Maxwell model moduli
clf;

f=linspace(0,10,1001);
E1=1*(f.^2)./(f.^2+1);
E2=(.5*f)./((.5*f).^2+1);

E_s=sqrt(E1.^2+E2.^2);

iE1=sqrt(E1.^2-E2.^2);
iE2=sqrt(E1.^2+E2.^2);
FontSizeNum=20;
iE_s=iE1/iE2;

figure(1)
subplot(1,2,1)
plot(f,0.96*E1,'linewidth',5,'Color',[0.4660 0.6740 0.1880]);
hold on
plot(f,E2,'--','linewidth',5,'Color',[0.9290 0.6940 0.1250]);
legend('$\frac{E_1}{E}$','$\frac{E_2}{E}$','Interpreter','latex','orientation','horizontal')
xlabel('Dimensionless Frequency ($\tau \omega$)','Interpreter','latex')
ylabel('Dimensionless Moduli')
grid on
ylim([0 1.4])
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)%,'fontWeight','bold')
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

subplot(1,2,2)
plot(f,E_s,'-.','linewidth',5,'Color',[0.6350 0.0780 0.1840]);
legend('$\frac{\sqrt{E_1^2+E_2^2}}{E}$','Interpreter','latex')
xlabel('Dimensionless Frequency ($\tau \omega$)','Interpreter','latex')
ylabel('Dimensionless Moduli')
grid on
ylim([0 1.4])
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)%,'fontWeight','bold')
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)

% plot of the moduli together
figure(2)
plot(f,0.96*E1,'linewidth',5,'Color',[0.4660 0.6740 0.1880]);
hold on
plot(f,E2,'--','linewidth',5,'Color',[0.9290 0.6940 0.1250]);
plot(f,E_s,'-.','linewidth',5,'Color',[0.6350 0.0780 0.1840]);
legend('$\frac{\textsf{E}_1}{\textsf{E}}$','$\frac{\textsf{E}_2}{\textsf{E}}$','$\frac{\sqrt{\textsf{E}_1^2+\textsf{E}_2^2}}{\textsf{E}}$','Interpreter','latex','orientation','horizontal')
xlabel('Dimensionless Frequency (\tau\omega)')
ylabel('Dimensionless Moduli')
grid on
ylim([0 1.4])
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',FontSizeNum)%,'fontWeight','bold')
set(gca,'FontName','Helvetica','FontSize',FontSizeNum)