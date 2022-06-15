clc; clear; close all; 

% Experiment parameters 
Nmonte            = 1;  % number of Monte carlo runs 
Tk                = 200; % number of samples 
x0                = [1; 10]; % inital state vector 
Ts                = 1;  % sampling time 
noise.process     = 1;  % process noise 
noise.measurement = 1;  % measurement noise  

% Space allocation for efficiency 
nx      = length(x0); 
nz      = 1; 
% System 
xk       = zeros(nx, Tk);     % true state vector 
zk       = zeros(1, Tk);      % measurement vector 
% Kalman filter
xk_hat   = zeros(nx, Tk);     % estiamted state 
Pk       = zeros(nx, nx, Tk); % estimated covariance 
% Error statistics 
error    = zeros(nx, Tk);    % error vector 
NIS      = zeros(1,Tk);      % NIS statistics 
NEES     = zeros(1, Tk);     % NEES statistics 
% Average error statistics 
avgNIS   = zeros(1,Tk); 
avgError = zeros(nx, Tk); 
avgNEES  = zeros(1,Tk); 


for n = 1:Nmonte
    for k = 1:Tk
        if k ==1 
            % Simulate system for k = 1 
            Model            = kinematic_models(Ts); 
            [xk(:, 1), F, Q] = Model.CWNA(x0, noise.process); 
            [zk(:,1), H, R]  = Model.position_measurements(xk(:,1), noise.measurement);  

            % Intialize stat_estimator class 
            Estimator        = state_estimator(F, Q, H, R); 

        else 
            % Simulate system for k = 2,...,Tk
            xk(:, k) = Model.CWNA(xk(:,k-1), noise.process);   
            zk(:, k) = Model.position_measurements(xk(:,k), noise.measurement);

            % Two point intialization using measurements from k = 1,2
            if k == 2 
                [xk_hat(:,1), Pk(:,:,1)]  = Estimator.two_point_differencing_WNA(zk(:,1:2), Ts);
                [error(:,k), NEES(:,k)] = Estimator.consistency_testing(xk(:,k), xk_hat(:,k), Pk(:,:,k)); 
            end
            
            % Kalman filter for k = 2,...,Tk 
            [xk_hat(:,k), Pk(:,:,k), NIS(:,k)] = Estimator.KF(xk_hat(:,k-1), Pk(:,:,k-1), zk(:,k)); 
              
            % Error statistics for k = 1,2...,Tk
            [error(:,1), NEES(:,1)] = Estimator.consistency_testing(xk(:,1), xk_hat(:,1), Pk(:,:,1));      
        end 
    end 

    avgError = avgError + error; 
    avgNIS   = avgNIS   + NIS; 
    avgNEES  = avgNEES  + NEES; 

end 


% Averaging for multiple-run tests ........................................
avgError = error/Nmonte; 
avgNEES  = NEES/Nmonte;
avgNIS   = NIS/Nmonte;

% single-run NIS ..........................................................
Lwindow  = 5; % window length
NIS      = movmean(NIS, Lwindow); % single-run time average NIS 

bounds.avgNEES = Estimator.confidence_bounds(nx, Nmonte);
bounds.avgNIS  = Estimator.confidence_bounds(nz, Nmonte); 
bounds.NIS     = Estimator.confidence_bounds(Lwindow, Nmonte); 

compliance.avgNEES = Estimator.calculate_compliance(avgNEES, bounds.avgNEES);
compliance.avgNIS  = Estimator.calculate_compliance(avgNIS, bounds.avgNIS);
compliance.NIS     = Estimator.calculate_compliance(NIS, bounds.NIS);



% plot parameters .........................................................



Plotter = plotter(Ts, Tk); 

Plotter.KF_estimate_plot(xk, xk_hat)

