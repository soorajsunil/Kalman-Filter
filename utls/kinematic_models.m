classdef kinematic_models

    % Discrete-time kinematic models for target tracking:
    %
    % Discretized continuous-time state space models
    %   1) CWNA
    %   2) CWPA
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
        sampling_period % sampling period in seconds
        model           % {CWNA, CWPA, DWNA, DWPA}
        process_noise
        measurement_noise
    end

    properties(Access=private)
        num_of_samples
    end

    methods

        function obj = kinematic_models(model, noise, sampling_period, num_of_samples)
            obj.model             = model;
            obj.sampling_period   = sampling_period;
            obj.num_of_samples    = num_of_samples;
            obj.process_noise     = noise.process;
            obj.measurement_noise = noise.measurement;
        end

        function [xk, zk, F, Q, H, R] = state_space_model(obj, x0)
            Ts  = obj.sampling_period;
            sigma_v = obj.process_noise;
            sigma_w = obj.measurement_noise;

            switch upper(obj.model)

                case 'CWNA' % Discretized continuous white noise acceleration (CWNA) model or second-order kinematic model (double integrator).
                    if ~iscolumn(x0) || size(x0,1) ~= 2
                        error('incorrect input state vector dimension (x0) try size 2x1')
                    end
                    F = [1 Ts;
                        0 1];
                    Q = [(Ts^3)/3 (Ts^2)/2;
                        (Ts^2)/2 Ts]*sigma_v;
                    H = [1 0];
                    Gamma = 1;

                case 'CWPA' % Discretized continuous Wiener process acceleration (CWPA) model or white noise jerk model.
                    if ~iscolumn(x0) || length(x0) ~= 3
                        error('incorrect input state vector dimension (x0) try size 3x1')
                    end
                    F = [1 Ts Ts^2/2;
                        0 1 Ts;
                        0 0 1];
                    Q = [(Ts^5)/20 (Ts^4)/8 (Ts^3)/6;
                        (Ts^4)/8 (Ts^3)/3 (Ts^2)/2;
                        (Ts^3)/6 (Ts^2)/2 Ts]*sigma_v;
                    H = [1 0 0];
                    Gamma = 1;

                case 'DWNA' % Discrete white noise acceleration (DWNA) model or piecewise constant white acceleration model
                    if ~iscolumn(x0) || length(x0) ~= 2
                        error('incorrect input state vector dimension (x0) try size 2x1')
                    end
                    F     = [1 Ts;
                        0 1];
                    Gamma = [Ts^2/2; Ts];
                    Q     = Gamma*(sigma_v^2)*Gamma';
                    H     = [1 0];

                case 'DWPA' % Discrete Wiener process acceleration (DWPA) model or piecewise constant Wiener process acceleration model
                    if ~iscolumn(x0) || length(x0) ~= 3
                        error('incorrect input state vector dimension (x0) try size 3x1')
                    end
                    F     = [1 Ts Ts^2/2;
                        0 1 Ts;
                        0 0 1];
                    Gamma = [Ts^2/2; Ts; 1];
                    Q     = Gamma*(sigma_v^2)*Gamma';
                    H     = [1 0 0];
                otherwise
                    error('select valid model name !!!')
            end
            R  = (sigma_w)^2; % measurement noise covariance
            xk = zeros(size(x0,1), obj.num_of_samples);
            zk = zeros(1, obj.num_of_samples);

            for k = 1:obj.num_of_samples
                if k == 1
                    [xk(:,1)] = x0 + Gamma*(sigma_v*randn);
                    [zk(:,1)] = linear_measurement_equation(xk(:,1), H, sigma_w);
                else
                    [xk(:,k)] = linear_state_equation(xk(:,k-1), F, sigma_v, Gamma);
                    [zk(:,k)] = linear_measurement_equation(xk(:,k), H, sigma_w);
                end
            end

        end
    end

end

function [xk] = linear_state_equation(xkm1, F, sigma_v, Gamma)
xk = F*xkm1 + Gamma*(sigma_v*randn); % process equation
end

function [zk] = linear_measurement_equation(xk, H, sigma_w)
zk = H*xk+ sigma_w*randn; % measurement equation
end


