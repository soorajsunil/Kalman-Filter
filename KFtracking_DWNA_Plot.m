close all; 

k = 1:system.samplingInterval:Tk; % x axis 

plt.SaveFigure   = 1;  
plt.Units        = 'centimeters'; 
plt.FontSize     = 10; 
plt.MarkerSize   = 6; 
plt.LineWidth    = 1.8; 
plt.Position     = [2, 8, 15, 10];
plt.FontName     = 'Arial';

%% ........................................................................
figure(Units=plt.Units, Position=plt.Position, Name='true vs estimate')

subplot(211) 
hold on; box on; grid on; 
plot(k, xk(1,:), 'k-', k, xk_hat(1,:), 'r--', 'Linewidth',plt.LineWidth) 
legend({'True','Estimated'}, 'location', 'northwest')
xlabel('Time (s)'); ylabel('Position (m)');
set(gca, 'fontsize', plt.FontSize,'FontName', plt.FontName)

subplot(212)
hold on; box on; grid on; 
plot(k, xk(2,:), 'k-', k, xk_hat(2,:), 'r--', 'Linewidth',plt.LineWidth) 
set(gca,'fontsize', plt.FontSize,'FontName',plt.FontName)
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend({'True','Estimated'}, 'location', 'northwest')

if plt.SaveFigure
    saveas(gca,[pwd '/DWNA_demo'],'epsc')
end

%% ........................................................................

figure(Units=plt.Units, Position=plt.Position,  Name='3 sigma bounds')

subplot(211)
y     = sqrt(squeeze(Pk(1,1,:))'); 
xbound = [k k(end:-1:1)];    
ybound = [-3.*y 3.*y(end:-1:1)];

bound           = fill(xbound,ybound, 'red');
bound.FaceColor = [1 0.8 0.8];      
bound.EdgeColor = 'none';  

hold on; box on; 
plot(k, avgError(1,:), 'r-', 'LineWidth', plt.LineWidth, 'MarkerSize', plt.MarkerSize)
set(gca, 'fontsize', plt.FontSize, 'FontName', plt.FontName)
xlabel('Time (s)')
ylabel('Position Error')
legend({'$ \pm 3 \sigma $ Estimated bounds', 'RMSE'}, 'interpreter', 'latex', Location='northeast')

subplot(212)
y              = sqrt(squeeze(Pk(2,2,:)))'; 
xbound          = [k k(end:-1:1)] ;         
ybound          = [-3.*y 3.*y(end:-1:1)];

bound           = fill(xbound,ybound,'red');
bound.FaceColor = [1 0.8 0.8];      
bound.EdgeColor = 'none';  

hold on; box on;  
plot(k, avgError(2,:)', 'r-', 'LineWidth', plt.LineWidth, 'MarkerSize', plt.MarkerSize)
legend({'$ \pm$ Estimated bounds', 'RMSE'}, 'interpreter', 'latex', Location='northeast')
xlabel('Time (s)'); ylabel('Velocity Error')
set(gca, 'fontsize', plt.FontSize, 'FontName', plt.FontName)

if plt.SaveFigure
    saveas(gca,[pwd '/DWNA_error_bounds'],'epsc')
end

%% ........................................................................

figure(Units=plt.Units, Position=plt.Position,  Name='average NIS')
xbound = [k k(end:-1:1)];    
ybound = [nis_r(1)*ones(1,Tk) nis_r(2)*ones(1,Tk)];

bound           = fill(xbound,ybound,'r');
bound.FaceColor = [1 0.8 0.8];         
bound.EdgeColor = 'none';  

hold on; box on; 
plot(k, avgNIS,'r-','Linewidth',plt.LineWidth)

set(gca, 'fontsize', plt.FontSize, 'FontName', plt.FontName)
xlabel('Time (s)');ylabel('Average NIS');
lgnd{2} = 'Average NIS'; lgnd{1} = [num2str(confidence) '% Probability region'];
legend(lgnd, 'location', 'northeast');
ylim([nis_r(1)-0.5 nis_r(2)+0.5])
xlim([1 Tk])
if plt.SaveFigure
    saveas(gca,[pwd '/DWNA_NIS'],'epsc')
end


%% ........................................................................

figure(Units=plt.Units, Position=plt.Position,  Name='average NEES')

xbound = [k k(end:-1:1)];    
ybound = [nees_r(1)*ones(1,Tk) nees_r(2)*ones(1,Tk)];

bound           = fill(xbound,ybound,'r');
bound.FaceColor = [1 0.8 0.8];         
bound.EdgeColor = 'none';  

hold on; box on; 
plot(k, avgNEES,'r-','Linewidth',plt.LineWidth)

set(gca, 'fontsize', plt.FontSize, 'FontName', plt.FontName)
xlabel('Time (s)');ylabel('Average NEES');
lgnd{2} = 'Average NEES'; lgnd{1} = [num2str(confidence) '% Probability region'];
legend(lgnd, 'location', 'northeast');
ylim([nees_r(1)-0.5 nees_r(2)+0.5])
xlim([1 Tk])
if plt.SaveFigure
    saveas(gca,[pwd '/DWNA_NEES'],'epsc')
end

%% ........................................................................

figure(Units=plt.Units, Position=plt.Position,  Name='average NMEE')

xbound = [k k(end:-1:1)];    
ybound = [nmee_r(1)*ones(1,Tk) nmee_r(2)*ones(1,Tk)];

subplot(211)
bound           = fill(xbound,ybound,'r');
bound.FaceColor = [1 0.8 0.8];         
bound.EdgeColor = 'none';  
hold on; box on; 
plot(k, avgNMEE(1,:),'r-','Linewidth',plt.LineWidth)
set(gca, 'fontsize', plt.FontSize, 'FontName', plt.FontName)
xlabel('Time (s)');ylabel('Average NMEE');
lgnd{2} = 'Average NMEE'; lgnd{1} = [num2str(confidence) '% Probability region'];
legend(lgnd, 'location', 'northeast');
ylim([nmee_r(1)-1 nmee_r(2)+1])
xlim([1 Tk])

subplot(212)
bound           = fill(xbound,ybound,'r');
bound.FaceColor = [1 0.8 0.8];         
bound.EdgeColor = 'none';  
hold on; box on; 
plot(k, avgNMEE(2,:),'r-','Linewidth',plt.LineWidth)
set(gca, 'fontsize', plt.FontSize, 'FontName', plt.FontName)
xlabel('Time (s)');ylabel('Average NMEE');
lgnd{2} = 'Average NMEE'; lgnd{1} = [num2str(confidence) '% Probability region'];
legend(lgnd, 'location', 'northeast');
ylim([nmee_r(1)-1 nmee_r(2)+1])
xlim([1 Tk])

if plt.SaveFigure
    saveas(gca,[pwd '/DWNA_NMEE'],'epsc')
end

clear lgnd xbound ybound
