classdef kinematic_models

    % Discrete-time linear kinematic models for target tracking:
    %   1) DWNA
    %   2) DWPA
    %
    % Code by: SOORAJ SUNIL
    % Date: 13 June 2022
    %
    % Reference:
    % [1]. Estimation with applications to tracking and navigation: theory
    % algorithms and software, Bar-Shalom, Yaakov and Li, X Rong and
    % Kirubarajan, Thiagalingam

    properties
        sampling_period % sampling period in seconds
        model           % {DWNA, DWPA}
        process_noise
        measurement_noise
        num_of_samples
    end

    methods
        function obj = kinematic_models(model, sampling_period, num_of_samples, process_noise, measurement_noise)
            obj.model             = model;
            obj.sampling_period   = sampling_period;
            obj.num_of_samples    = num_of_samples;
            obj.process_noise     = process_noise;
            obj.measurement_noise = measurement_noise;
        end
        function [xk, zk, F, Q, H, R] = simulate(obj, x0)
            Ts  = obj.sampling_period;
            sigma_v = obj.process_noise;
            sigma_w = obj.measurement_noise;
            switch upper(obj.model)
                case 'DWNA'
                    % Discrete white noise acceleration (DWNA) model or piecewise constant white acceleration model
                    if ~iscolumn(x0) || length(x0) ~= 2
                        error('incorrect input state vector dimension (x0) try size 2x1')
                    end
                    F = [1 Ts;
                        0 1];                     % state transition matrix
                    Gamma = [Ts^2/2; Ts];         % noise gain matrix
                    Q = Gamma*(sigma_v^2)*Gamma'; % noise covariance
                    H = [1 0];                    % measurement matrix
                case 'DWPA'
                    % Discrete Wiener process acceleration (DWPA) model or piecewise constant Wiener process acceleration model
                    if ~iscolumn(x0) || length(x0) ~= 3
                        error('incorrect input state vector dimension (x0) try size 3x1')
                    end
                    F = [1 Ts Ts^2/2;
                        0 1 Ts;
                        0 0 1];                    % state transition matrix
                    Gamma = [Ts^2/2; Ts; 1];       % noise gain matrix
                    Q = Gamma*(sigma_v^2)*Gamma';  % noise covariance
                    H = [1 0 0];                   % measurement matrix
                otherwise
                    error('invalid model name!')
            end
            R  = (sigma_w)^2;                           % measurement noise covariance
            xk = zeros(size(x0,1), obj.num_of_samples); % state vector 
            zk = zeros(1, obj.num_of_samples);          % measurement vector
            for k = 1:obj.num_of_samples
                if k == 1
                    [xk(:,1)] = x0 + Gamma*(sigma_v*randn);
                    [zk(:,1)] = obj.measurement_eq(xk(:,1), H, sigma_w);
                else
                    [xk(:,k)] = obj.process_eq(xk(:,k-1), F, Gamma, sigma_v);
                    [zk(:,k)] = obj.measurement_eq(xk(:,k), H, sigma_w);
                end
            end
        end
    end
    methods(Static)
        function [xk] = process_eq(xkm1, F, Gamma, sigma_v)
            xk = F*xkm1 + Gamma*(sigma_v*randn); % process equation
        end
        function [zk] = measurement_eq(xk, H, sigma_w)
            zk = H*xk+ sigma_w*randn; % measurement equation
        end
    end
end



