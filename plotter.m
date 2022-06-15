classdef plotter
    
    properties
        SaveFigure   = 0;  
        Units        = 'centimeters'; 
        FontSize     = 12; 
        MarkerSize   = 8; 
        LineWidth    = 1.8; 
        Position     = [2, 8, 15, 10];
        FontName     = 'Arial';
        k            

    end
    
    methods

        function obj = plotter(Ts, Tk)
            obj.k = (1:1:Tk).*Ts; % x axis (time); 
        end 
        
        function KF_estimate_plot(obj, xk, xk_hat)
            figure(Units=obj.Units, Position=obj.Position, Name='True vs estimate')
            [nx, ~] = size(xk);
            for p = 1:nx 
                subplot(nx, 1, p) 
                hold on; box on; grid on; 
                set(gca, 'fontsize', obj.FontSize,'FontName', obj.FontName)
                plot(obj.k, xk(p,:), 'k-', obj.k, xk_hat(p,:), 'r--', 'Linewidth',obj.LineWidth) 
                xlabel('Time (s)'); 
                legend({'True','Estimated'}, 'location', 'northwest')

                switch p 
                    case 1 
                        ylabel('Position (m)');
                    case 2 
                        ylabel('Velocity (m/s)')
                    case 3 
                        ylabel('Acceleration (m/s^2)')
                end
            end 
            if obj.SaveFigure
                saveas(gca,[pwd '/KF_estimate'],'epsc')
            end
        end
    end
end

