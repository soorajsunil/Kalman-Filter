% time vector 
time_k = (0:Nsamples-1)*Ts; % x axis 
time_k = time_k(3:end); 

% confidence interval 
Nz = size(zk,1);
[nis_bounds, nees_bounds, ~, confidence] = confidence_bounds(Nx, Nz, Nmonte); 

% plot properties 
plt.save         = true; 
plt.saveFormat   = 'bmp'; 
plt.Units        = 'centimeters'; 
plt.FontSize     = 10; 
plt.MarkerSize   = 6; 
plt.LineWidth    = 1.8; 
plt.Position     = [2, 8, 15, 10];
plt.FontName     = 'Arial';

figure(Units=plt.Units, Position=plt.Position, Name='true vs estimate')
for kx = 1:Nx
    subplot(Nx,1,kx)
    hold on; box on; grid on;
    plot(time_k, xk(kx,:), 'k-', time_k, xk_hat(kx,:), 'r--', 'Linewidth',plt.LineWidth)
    legend({'True','Estimated'}, 'location', 'northwest')
    xlabel('Time (s)'); 
    set(gca, 'fontsize', plt.FontSize,'FontName', plt.FontName)
    axis('tight')

    switch kx 
        case 1
            ylabel('Position (m)');

        case 2 
            ylabel('Velocity (m/s)')
        case 3 
            ylabel('Acceleration (m/s^2)')
    end 
end 
if plt.save
    saveas(gca,fullfile('plots/', strcat(model,'state')),plt.saveFormat)
end 

figure(Units=plt.Units, Position=plt.Position,  Name='3 sigma bounds')
for kx = 1:Nx
    subplot(Nx,1,kx)
    y     = sqrt(squeeze(Pk(kx,kx,:))');
    xbound = [time_k time_k(end:-1:1)];
    ybound = [-3.*y 3.*y(end:-1:1)];
    bound           = fill(xbound,ybound, 'red');
    bound.FaceColor = [1 0.8 0.8];
    bound.EdgeColor = 'none';
    hold on; box on;
    plot(time_k, avgError(kx,:), 'r-', 'LineWidth', plt.LineWidth, 'MarkerSize', plt.MarkerSize)
    set(gca, 'fontsize', plt.FontSize, 'FontName', plt.FontName)
    xlabel('Time (s)')
    ylabel('Position Error')
    legend({'$ \pm 3 \sigma $ Estimated bounds', 'RMSE'}, 'interpreter', 'latex', Location='northeast')
    axis('tight')
     switch kx 
        case 1
            ylabel('Position (m)');

        case 2 
            ylabel('Velocity (m/s)')
        case 3 
            ylabel('Acceleration (m/s^2)')
    end 
end 
if plt.save
    saveas(gca,fullfile('plots/', strcat(model,'covariance')),plt.saveFormat)
end 

figure(Units=plt.Units, Position=plt.Position,  Name='average NIS')
xbound = [time_k time_k(end:-1:1)];    
ybound = [nis_bounds(1)*ones(1,K) nis_bounds(2)*ones(1,K)];
bound           = fill(xbound,ybound,'r');
bound.FaceColor = [1 0.8 0.8];         
bound.EdgeColor = 'none';  
hold on; box on; 
plot(time_k, avgNIS,'r-','Linewidth',plt.LineWidth)
set(gca, 'fontsize', plt.FontSize, 'FontName', plt.FontName)
xlabel('Time (s)');ylabel('Average NIS');
lgnd{2} = 'Average NIS'; lgnd{1} = [num2str(confidence) '% Probability region'];
legend(lgnd, 'location', 'northeast');
ylim([nis_bounds(1)-0.5 nis_bounds(2)+0.5])
xlim([1 Nsamples])
if plt.save
    saveas(gca,fullfile('plots/', strcat(model,'nis')),plt.saveFormat)
end 

figure(Units=plt.Units, Position=plt.Position,  Name='average NEES')
xbound = [time_k time_k(end:-1:1)];    
ybound = [nees_bounds(1)*ones(1,K) nees_bounds(2)*ones(1,K)];
bound           = fill(xbound,ybound,'r');
bound.FaceColor = [1 0.8 0.8];         
bound.EdgeColor = 'none';  
hold on; box on; 
plot(time_k, avgNEES,'r-','Linewidth',plt.LineWidth)
set(gca, 'fontsize', plt.FontSize, 'FontName', plt.FontName)
xlabel('Time (s)');ylabel('Average NEES');
lgnd{2} = 'Average NEES'; lgnd{1} = [num2str(confidence) '% Probability region'];
legend(lgnd, 'location', 'northeast');
ylim([nees_bounds(1)-0.5 nees_bounds(2)+0.5])
xlim([1 Nsamples])
if plt.save
    saveas(gca,fullfile('plots/', strcat(model,'nees')),plt.saveFormat)
end 

clear lgnd xbound ybound kx bound y
