classdef state_estimator
    % Linear Kalman filter 

    % Reference:
    % [1]. Estimation with applications to tracking and navigation: theory 
    % algorithms and software, Bar-Shalom, Yaakov and Li, X Rong and 
    % Kirubarajan, Thiagalingam
 
    properties
        F % state transition matrix 
        Q % process noise covariance 
        H % measurement matrix  
        R % measurement noise covariance 
    end
    
    methods
        function obj = state_estimator(F, Q, H, R)
          obj.F = F;
          obj.Q = Q; 
          obj.H = H;
          obj.R = R; 
        end

        function [xk0_hat, P0] = TPD_initialization_WNA(obj, zk, Ts)
            % Two point differencing filter initialization for second-order
            % models [Section 5.5.3]
            xk0_hat = [zk(:,1); (zk(:,2)-zk(:,1))/Ts];         
            P0      = [obj.R -obj.R/Ts; -obj.R/Ts 2*obj.R/Ts^2]; 
        end 
    
        function [xk_hat, Pk, NIS] = Kalman_filter(obj, xkm1_hat, Pkm1, zk)
          
            [nx, ~] = size(xkm1_hat); 
            [~, Tk] = size(zk); 
            xk_hat = zeros(nx, Tk); 
            Pk     = zeros(nx, nx, Tk); 
            NIS    = zeros(1,Tk); 

            if Tk == 1 
                [xk_hat, Pk, NIS] = KF(obj, xkm1_hat, Pkm1, zk);      
            else
                xk_hat(:,1) = xkm1_hat;
                Pk(:,:,1)   = Pkm1; 
                for k = 2:Tk 
                    [xk_hat(:,k), Pk(:,:,k), NIS(:,k)] = obj.KF_cycle(xk_hat(:,k-1), Pk(:,:,k-1), zk(:,k));      
                end
            end 

        end

        function [xk_hat, Pk, NIS] = KF_cycle(obj, xkm1_hat, Pkm1, zk)
            % Kalman filter one cycle [Section 5.2.4]

            % State and covariance prediction 
            xpred  = obj.F*xkm1_hat;            % predicted state 
            Ppred  = obj.F*Pkm1*obj.F' + obj.Q; % predicted covariance                 
            % Measurement prediction            
            zk_hat = obj.H*xpred;                % predicted measurement 
            inov   = zk - zk_hat;                % residual 
            S      = obj.R + obj.H*Ppred*obj.H'; % innovation covariance
            NIS    = inov'*S^(-1)*inov;   % normalized innovation squared
            % State and covariance update 
            G      = Ppred*obj.H'/(S); % kalman gain
            xk_hat = xpred + G*inov;   % updated state
            % Pk   = Ppred - G*S*G';   % updated covariance   
            % Joseph form (alternative covarian ce update for numerical stability)
            Pk     = (eye(size(obj.Q,1)) - G*obj.H)*Ppred*(eye(size(obj.Q,1)) ...
                        - G*obj.H)' + G*obj.R*G'; 
        end

        function [xk_error] = estimation_error(~, xk, xk_hat)
            % estimation error
            xk_error = xk - xk_hat; 
        end 

        function [NEES] = NEES(~, xk_error, Pk)
            % normalized estimation error squared
             
            [~, Tk] = size(xk_error); 
            NEES    = zeros(1, Tk); 
 
            for k = 1:Tk
                NEES(:,k)  = xk_error(:,k)'*Pk(:,:,k)^(-1)*xk_error(:,k); 
            end
        end 
    end
end

