% Tracking the dynamic state (position and velocity) of a
% constant velocity target using the Kalman filter.
% Writer: SOORAJ SUNIL
% Date: 10 May 2022
% References:
% [1]. Estimation with applications to tracking and navigation: Theory
% algorithms and software, Bar-Shalom, Yaakov and Li, X Rong and
% Kirubarajan, Thiagalingam

% START OF MAIN ... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; 
close all;

% Experiment parameters
Nsamples = 100;       % number of samples
dt       = 0.5;       % sampling period
x0       = [420; 5];  % initial true state
sigma_v  = 1;         % process noise std. deviation
sigma_w  = 1;         % measurement noise std. deviation

% Simulate system (* additional two samples to initialize the filter)
[xk, zk, F, H, Q, R] = DWNAmodel(x0,dt,(Nsamples+2),sigma_v,sigma_w);

% Kalman filter
Nx    = length(x0);
xkhat = zeros(Nx, Nsamples);      % estimated state vector
Pk    = zeros(Nx, Nx, Nsamples);  % estimated covariance
NIS   = zeros(1, Nsamples);       % normalized innovation squared (NIS)

% Filter initialization
[x0hat, P0, xk, zk] = TPDintialization(xk, zk, R, dt);

for k = 1:Nsamples
    [xkhat(:,k), Pk(:,:,k), NIS(:,k)] = KalmanFilter(x0hat, P0, F, H, Q, R, zk(:,k)); % Kalman filter
    x0hat = xkhat(:,k); % update previous estimate
    P0    = Pk(:,:,k);  % update previous covariance
end
clear k P0 x0hat

% Figures:
figure(Name='Position'); hold on;
k = (1:1:Nsamples).*dt; % plot x-axis
plot(k, xk(1,:), 'k-', k, xkhat(1,:), 'r--', 'Linewidth',2)
xlabel('Time (s)'); ylabel('Position (m)');
legend({'True','KF Estimate'}); hold off;

figure(Name='Velocity'); hold on;
plot(k, xk(2,:), 'k-', k, xkhat(2,:), 'r--', 'Linewidth',2)
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend({'True','KF Estimate'}); hold off;

% END OF MAIN ... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x0_hat,P0,xk,zk] = TPDintialization(xk,zk,R,T)
% two-point differencing using the first two measurements
x0_hat = [zk(1); (zk(2)-zk(1))/T];
P0     = [1 -1/T; -1/T 2/T^2].*(R)^2;
% remove the two points from the true system and measurement vectors
xk = xk(:,3:length(xk));
zk = zk(:,3:length(zk));
end

function [xk, zk, F, H, Q, R] = DWNAmodel(x0,dt,Tk,sigma_v,sigma_w)
% calculate system matrices
F     = [1 dt; 0 1];
Gamma = [0.5*(dt^2); dt];
H     = [1 0];
Q     = Gamma*((sigma_v)^2)*Gamma';
R     = (sigma_w)^2;
% Initialize:
xk(:,1) = x0 + (sigma_v).*randn(size(x0));
zk(:,1) = H*xk(:,1) + sigma_w*randn;
for k = 2:Tk
    xk(:,k) = F*xk(:,k-1) + Gamma*(sigma_v*randn); % process equation
    zk(:,k) = H*xk(:,k) + sigma_w*randn;           % measurement equation
end
end

function [xkhat, Pk, NIS] = KalmanFilter(x0, P0, F, H, Q, R, zk)
% State and covariance prediction
xpred  = F*x0;                % predicted state
Ppred  = F*P0*F' + Q;         % predicted covariance
% Measurement prediction
zk_hat = H*xpred;             % predicted measurement
inov   = zk - zk_hat;         % residual
S      = R + H*Ppred*H';      % innovation covariance
NIS    = inov'/S*inov;   % normalized innovation squared
% State and covariance update
G      = Ppred*H'/(S);        % kalman gain
xkhat = xpred + G*inov;       % updated state
% Pk     = Ppred - G*S*G';    % updated covariance
% Joseph form (alternative covariance update for numerical stability)
Pk     = (eye(size(Q,1)) - G*H)*Ppred*(eye(size(Q,1)) - G*H)' + G*R*G';
end
