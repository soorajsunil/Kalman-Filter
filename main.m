clc; clear; close all; 

% Experiment parameters: 
Nmonte = 1;   % number of Monte carlo runs 
Tk     = 40; % number of samples 

% System parameters:
x0                = [1; 10]; % inital state vector 
Ts                = 1;       % sampling time 
noise.process     = 1;       % process noise 
noise.measurement = 1;       % measurement noise  

% Space allocation for efficiency:  
nx     = length(x0);     % state vector dimension 
nz     = 1;              % measurement vector dimension 
xk     = zeros(nx, Tk);  % true state vector 
zk     = zeros(1, Tk);   % measurement vector 




for n = 1:Nmonte

    for k = 1:Tk

        if k == 1 
            % Simulate system for k = 1 
            Model           = kinematic_models(Ts); 
            [xk(:,1), F, Q] = Model.CWNA(x0, noise.process); 
            [zk(:,1), H, R] = Model.position_measurements(xk(:,1), noise.measurement);  
            Estimator       = state_estimator(F, Q, H, R); 

        else 
            % Simulate system for k = 2,...,Tk
            xk(:,k) = Model.CWNA(xk(:,k-1), noise.process);   
            zk(:,k) = Model.position_measurements(xk(:,k), noise.measurement);

        end 

    end 
    clear k 
    % Intiialize filter using two measurement samples 
    [x0_hat, P0] = Estimator.TPD_initialization_WNA(zk(:,1:2), Ts); 
    
    % Kalman filter (batch)
    [xk_hat, Pk, NIS] = Estimator.Kalman_filter(x0_hat, P0, zk);

    % Estimation error 
    [xk_error]   = Estimator.estimation_error(xk, xk_hat);

    % NEES statistics 
    [NEES]       = Estimator.NEES(xk_error, Pk); 

end 
clear n 


plot(xk_error(1,:))

