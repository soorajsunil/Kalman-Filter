classdef kinematic_models
    % Discrete-time kinematic models: 
    %
    % Discretized continuous-time state space models
    %   1) discretized CWNA
    %   2) discretized CWPA
    % Discrete-time state space models
    %   3) DWNA
    %   4) DWPA 
    %
    % Code by : SOORAJ SUNIL
    % https://scholar.google.ca/user=mX2uugYAAAAJ&hl=en
    % Date    : 13 June 2022
    %
    % Reference:
    % [1]. Estimation with applications to tracking and navigation: theory 
    % algorithms and software, Bar-Shalom, Yaakov and Li, X Rong and 
    % Kirubarajan, Thiagalingam
  
    properties
        samplingTime 
    end

    methods
        function obj = kinematic_models(samplingTime)
            obj.samplingTime  = samplingTime; 
        end 

        function [xk, F, Q] = CWNA(obj, xkm1, qtild)
            % Discretized continuous white noise acceleration (CWNA) model 
            % or second-order kinematic model (double integrator).
            if ~iscolumn(xkm1) || length(xkm1) ~= 2 
                error(['incorrect state vector dimension;' ...
                    ' try column vector of size 2x1'])
            end 
            
            T  = obj.samplingTime; 
            F  = [1 T; 0 1];                       % state transition matrix 
            Q  = [T^3/3 T^2/2; T^2/2 T]*qtild;     % noise covariance matrix          
            xk = F*xkm1 + qtild*randn(size(xkm1)); % state equation 
        end

        function [xk, F, Q] = CWPA(obj, xkm1, qtild)
            % Discretized continuous Wiener process acceleration (CWPA) model
            % or white noise jerk model.
            if ~iscolumn(xkm1) || length(xkm1) ~= 3 
                error(['incorrect state vector dimension;' ...
                    ' try column vector of size 3x1'])
            end 

            T  = obj.samplingTime; 
            F  = [1 T T^2/2; 0 1 T; 0 0 1];        % state transition matrix 
            Q  = [T^5/20 T^4/8 T^3/6;
                  T^4/8 T^3/3 T^2/2; 
                  T^3/6 T^2/2 T]*qtild;            % noise covariance matrix   
    
            xk = F*xkm1 + qtild*randn(size(xkm1)); % state equation 
        end

        function [xk, F, Q] = DWNA(obj, xkm1, sigma_v)
            % Discrete white noise acceleration (DWNA) model
            % or piecewise constant white acceleration model
            if ~iscolumn(xkm1) || length(xkm1) ~= 2 
                error(['incorrect state vector dimension;' ...
                    ' try column vector of size 2x1'])
            end 

            T     = obj.samplingTime; 
            F     = [1 T; 0 1];               % state transition matrix 
            Gamma = [T^2/2; T];               % noise gain matrix
            Q     = Gamma*(sigma_v^2)*Gamma'; % noise covariance matrix          
            xk    = F*xkm1 + Gamma*(sigma_v*randn); % state equation 
        end

        function [xk, F, Q] = DWPA(obj, xkm1, sigma_v)
            % Discrete Wiener process acceleration (DWPA) model,
            % or piecewise constant Wiener process acceleration model 
            if ~iscolumn(xkm1) || length(xkm1) ~= 3 
                error(['incorrect state vector dimension;' ...
                    ' try column vector of size 3x1'])
            end 

            T     = obj.samplingTime; 
            F     = [1 T T^2/2; 0 1 T; 0 0 1];      % state transition matrix 
            Gamma = [T^2/2; T; 1];                  % noise gain matrix
            Q     = Gamma*(sigma_v^2)*Gamma';       % noise covariance matrix          
            xk    = F*xkm1 + Gamma*(sigma_v*randn); % state equation 
        end

        function [zk, H, R] = position_measurements(~, xk, sigma_w)
            % generate position measurements        
            if ~iscolumn(xk)
                error('state must be a column vector')
            end 
            switch length(xk) 
                case 2 
                    H = [1 0]; % second-order (CWNA and DWNA models) 
                case 3 
                    H = [1 0 0]; % third-order (CWPA and DWPA models)  
                otherwise 
                    error('state must be 2x1 or 3x1 vector')
            end 
            R  = sigma_w^2; % measurement noise covariance 
            zk = H*xk + sigma_w*randn; % measurement equation 
        end 
    end
end

