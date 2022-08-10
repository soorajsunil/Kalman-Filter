% % Tracking the dynamic state (position, velocity and acceleration) using the Kalman filter
% Code by : SOORAJ SUNIL
% Date    : 10 May 2022
%
% Reference:
% [1]. Estimation with applications to tracking and navigation: theory
% algorithms and software, Bar-Shalom, Yaakov and Li, X Rong and
% Kirubarajan, Thiagalingam
clc; clear; close all; addpath('utls/')

% % Experiment parameters
model    = 'DWNA';
Nsamples = 100;  % experiment time duration
Nmonte   = 50;   % number of monte carlo runs

% % System parameters
Ts = 1;
x0 = [0; 10];             % mean of intial state
noise.process     = 1;    % process noise std. deviation
noise.measurement = 1;    % measurement noise std. deviation
system = kinematic_models(model, noise, Ts, Nsamples); % intialize class

for n = 1:Nmonte
    % % Simulate system 
    [xk, zk, F, Q, H, R] = system.state_space_model(x0);

    % % Filter initialization (note: two samples are removed to initialize the filter)
    [xk, zk, x0_hat, P0] = two_point_differencing(model, xk, zk, Ts, R);

    % % Kalman filter
    [xk_hat, Pk, nis]    = kalman_filter(x0_hat, P0,zk, F, Q, H, R);

    % % Monte carlo averaging
    if Nmonte > 1
        if n == 1 % allocate space
            [Nx, K] = size(xk);
            avgError = zeros(Nx, K);  % average estimation error
            avgNIS   = zeros(1, K);   % average NIS
            avgNEES  = zeros(1, K);   % average normalized estimation error squared
            avgNMEE  = zeros(Nx, K);  % average mean estimation error
        end
        fprintf('Monte carlo no. = %d / %d \n ',n, Nmonte);
        % Error statistics
        [xk_error, NEES, NMEE] = error_statistics(xk, xk_hat, Pk);
        avgError = avgError + xk_error;
        avgNIS   = avgNIS   + nis;
        avgNEES  = avgNEES  + NEES;
        avgNMEE  = avgNMEE  + NMEE;
        clear NMEE NEES bias squaredError

        if n == Nmonte % compute average
            avgError = avgError/Nmonte;
            avgNIS   = avgNIS/Nmonte;
            avgNEES  = avgNEES/Nmonte;
            avgNMEE  = avgNMEE/Nmonte;
            break; 
        end
    end
end 
clear n
