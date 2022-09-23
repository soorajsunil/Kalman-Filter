clc; clear; close all; addpath('utls/')

% % Experiment parameters
model    = 'DWNA';
Nsamples = 100;  % experiment time duration
Nmonte   = 50;   % number of monte carlo runs

% % System parameters
Ts = 1;                   % sampling period
x0 = [0; 10];             % mean of intial state
process_noise     = 1;    % process noise std. deviation
measurement_noise = 1;    % measurement noise std. deviation

% intialize class
system = kinematic_models(model, Ts, Nsamples, process_noise, measurement_noise);

for n = 1:Nmonte
    % Simulate system
    [xk, zk, F, Q, H, R] = system.simulate(x0);

    % Filter initialization (note: two samples are removed to initialize the filter)
    [xk, zk, x0_hat, P0] = two_point_differencing(model, xk, zk, Ts, R);

    % Kalman filter
    [xk_hat, Pk, nis] = kalman_filter(x0_hat, P0, zk, F, Q, H, R);

    % Monte carlo averaging
    if Nmonte > 1
        if n == 1 % allocate space
            [Nx, K]   = size(xk);
            avg_error = zeros(Nx,K);  % average estimation error
            avg_NIS   = zeros(1,K);   % average NIS
            avg_NEES  = zeros(1,K);   % average normalized estimation error squared
            avg_NMEE  = zeros(Nx,K);  % average mean estimation error
        end
        fprintf('Monte carlo no. = %d / %d \n ',n, Nmonte);
        % Error statistics
        [xk_error, NEES, NMEE] = error_statistics(xk, xk_hat, Pk);
        avg_error = avg_error + xk_error;
        avg_NIS   = avg_NIS   + nis;
        avg_NEES  = avg_NEES  + NEES;
        avg_NMEE  = avg_NMEE  + NMEE;
        clear NMEE NEES bias squaredError
        if n == Nmonte % compute average
            avg_error = avg_error/Nmonte;
            avg_NIS   = avg_NIS/Nmonte;
            avg_NEES  = avg_NEES/Nmonte;
            avg_NMEE  = avg_NMEE/Nmonte;
            break;
        end
    end
end
clear n

%% Plot 
plot_KFdemo
