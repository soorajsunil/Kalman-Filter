classdef state_estimator
 
    properties
        F % state transition matrix 
        Q % process noise covariance 
        H % measurement matrix  
        R % measure noise covariance 
    end
    
    methods
        function obj = state_estimator(F, Q, H, R)
          obj.F = F;
          obj.Q = Q; 
          obj.H = H;
          obj.R = R; 
        end

        function [xk0_hat, P0] = two_point_differencing_WNA(obj, zk, Ts)
            % Filter intialization 
            xk0_hat = [zk(:,1); (zk(:,2)-zk(:,1))/Ts];           % Eq. 5.5.3-4
            P0      = [obj.R -obj.R/Ts; -obj.R/Ts 2*obj.R/Ts^2]; % Eq. 5.5.3-5
        end 
    
        function [xk_hat, Pk, nis] = KF(obj, xkm1_hat, Pkm1, zk)
            % Kalman filter  
            % State and covariance prediction 
            xpred  = obj.F*xkm1_hat;            % predicted state 
            Ppred  = obj.F*Pkm1*obj.F' + obj.Q; % predicted covariance                 
            % Measurement prediction            
            zk_hat = obj.H*xpred;                % predicted measurement 
            inov   = zk - zk_hat;                % residual 
            S      = obj.R + obj.H*Ppred*obj.H'; % innovation covariance
            nis    = inov'*S^(-1)*inov;   % normalized innovation squared
            % State and covariance update 
            G      = Ppred*obj.H'/(S); % kalman gain
            xk_hat = xpred + G*inov;   % updated state
            % Pk   = Ppred - G*S*G';   % updated covariance   
            % Joseph form (alternative covariance update for numerical stability)
            Pk     = (eye(size(obj.Q,1)) - G*obj.H)*Ppred*(eye(size(obj.Q,1)) ...
                        - G*obj.H)' + G*obj.R*G'; 

        end

        function [xk_tild, nees] = consistency_testing(~, xk, xk_hat, Pk)
            xk_tild = xk - xk_hat;              % estimation error
            nees    = xk_tild'*Pk^(-1)*xk_tild; % normalized estimation error squared 
        end 

        function [bounds, confidence] = confidence_bounds(~, dof, Nmonte)   
            alpha   = 0.05;
            b1      = chi2inv(alpha/2, Nmonte*dof) /Nmonte;
            b2      = chi2inv(1-(alpha/2), Nmonte*dof) /Nmonte;
            bounds  = [b1; b2]; 
            confidence = 100*(1-alpha);
        end 

        function [compliance] = calculate_compliance(~, test, bounds)
            Noutliers  = length(find(test>bounds(2))) + length(find(test<bounds(1))); 
            compliance = 100 - 100*Noutliers/length(test);
            compliance = abs(compliance); 
        end 

    end
end

