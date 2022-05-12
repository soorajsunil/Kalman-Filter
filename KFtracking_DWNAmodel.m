clc; clear; close all; 

% Tracking the dynamic state (position and velocity) of a
% discrete white noise acceleration (DWNA) model using the Kalman filter 

% Code by : SOORAJ SUNIL
% https://scholar.google.ca/user=mX2uugYAAAAJ&hl=en
% Date    : 10 May 2022 
% Reference:
% [1]. Estimation with applications to tracking and navigation: theory 
% algorithms and software, Bar-Shalom, Yaakov and Li, X Rong and 
% Kirubarajan, Thiagalingam

% Experiment parameters ...................................................

Tk     = 100;   % experiment time duration 
Nmonte = 50;   % number of monte carlo runs 

% System parameters .......................................................
system.samplingInterval = 1;              
system.intialState      = [0; 10];  % mean of intial state        
noise.process           = 1;        % process noise std. deviation 
noise.measurement       = 1;        % measurement noise std. deviation 

% .........................................................................
    
% Space allocation 
Nx       = length(system.intialState);   
avgError = zeros(Nx, Tk);      % average estimation error 
avgNIS   = zeros(1, Tk);       % average NIS
avgNEES  = zeros(1,Tk);        % average normalized estimation error squared
avgNMEE  = zeros(Nx,Tk);       % average mean estimation error 

% .........................................................................

for n = 1:Nmonte
    fprintf('Monte carlo no. = %d / %d \n ',n, Nmonte);

    % Simulate system (* additional two samples to initialize the filter) 
    [xk, zk, F, H, Q, R] = simulate_DWNA_system(noise, system, (Tk+2));

    % Filter initialization
    [xk, zk, xkm1_hat, Pkm1] = ...
        two_point_differencing(xk, zk, R, system.samplingInterval);

    % Kalman filter 
    xk_hat  = zeros(Nx, Tk);      % estimated state vector
    Pk      = zeros(Nx, Nx, Tk);  % estimated covariance
    NIS     = zeros(1, Tk);       % normalized innovation squared (NIS) 

    for k = 1:Tk
        [xk_hat(:,k), Pk(:,:,k), NIS(:,k)] = ...
            kalman_filter(xkm1_hat, Pkm1, F, H, Q, R, zk(:,k)); 
        xkm1_hat = xk_hat(:,k);  % update previous estimate
        Pkm1     = Pk(:,:,k);    % update previous covariance 
    end 
    clear  k pkm1 xkm1_hat
    
    % Error statistics 
    [xk_tild, NEES, NMEE] = error_statistics(xk, xk_hat, Pk); 

    avgError = avgError + xk_tild; 
    avgNIS   = avgNIS   + NIS; 
    avgNEES  = avgNEES  + NEES; 
    avgNMEE  = avgNMEE  + NMEE; 
    clear NMEE NEES bias squaredError 

end 
clear n 

% calculate average ........................................... ............
avgError = avgError./Nmonte; 
avgNEES  = avgNEES./Nmonte;
avgNMEE  = avgNMEE./Nmonte; 
avgNIS   = avgNIS./Nmonte;

% confidence interval 
[Nz, ~] = size(zk);
[nis_r, nees_r, nmee_r, confidence] = confidence_interval(Nx, Nz, Nmonte); 


% Plot ....................................................................

KFtracking_DWNA_Plot

% .........................................................................

function [xk, zk, x0_hat, P0] = two_point_differencing(xk, zk, R, T)

% two point differencing using first two measurements 
x0_hat = [zk(1); (zk(2)-zk(1))/T]; 
P0     = [1 -1/T; 
          -1/T 2/T^2].*(R)^2; 

% remove the two-points from the true sytsem and measurement vectors 
xk     = xk(:, 3:length(xk)); 
zk     = zk(:, 3:length(zk)); 

end 

function [xk, zk, F, H, Q, R] = simulate_DWNA_system(noise, system, Tk)

% calculate system matrices
F      = [1 system.samplingInterval; 
          0 1];      
Gamma  = [0.5*(system.samplingInterval^2); 
          system.samplingInterval]; 
H      = [1 0]; 

% calculate noise co/variances 
Q       = Gamma*(noise.process^2)*Gamma'; 
R       = (noise.measurement)^2;     

xk(:,1) = system.intialState + (noise.process).*randn(size(system.intialState)); 
zk(:,1) = H*xk(:,1) + noise.measurement*randn; 

for k = 2:Tk
    xk(:,k) = F*xk(:,k-1) + Gamma*(noise.process*randn); % process equation
    zk(:,k) = H*xk(:,k) + noise.measurement*randn; % measurement equation 
end

end 

function [xk_hat, Pk, nis] = kalman_filter(xkm1_hat, Pkm1, F, H, Q, R, zk)

    % State and covariance prediction 
    xpred  = F*xkm1_hat;    % predicted state 
    Ppred  = F*Pkm1*F' + Q; % predicted covariance
            
    % Measurement prediction            
    zk_hat = H*xpred;             % predicted measurement 
    inov   = zk - zk_hat;         % residual 
    S      = R + H*Ppred*H';      % innovation covariance
    nis    = inov'*inv(S)*inov; % normalized innovation squared

    % State and covariance update 
    G      = Ppred*H'/(S);     % kalman gain
    xk_hat = xpred + G*inov;   % updated state
%   Pk     = Ppred - G*S*G';   % updated covariance   

    % Joseph form (alternative covariance update for numerical stability)
    Pk     = (eye(size(Q,1)) - G*H)*Ppred*(eye(size(Q,1)) - G*H)' + G*R*G'; 
end 

function [xk_tild, NEES, NMEE] = error_statistics(xk, xk_hat, Pk)

[Nx, Tk] = size(xk);

% estimation error
xk_tild = xk - xk_hat; 

% normalized estimation error squared 
NEES = zeros(1,Tk); 
for k = 1:Tk
    NEES(k) = xk_tild(:,k)'*(Pk(:,:,k))^(-1)*xk_tild(:,k); 
end 

% normalized mean estimation error
NMEE   = zeros(Nx,Tk); 
for k = 1:Tk
    for i = 1:Nx
       NMEE(i,k) = (xk_tild(i,k)/sqrt(Pk(i,i,k)));
    end
end
      
end 

function [nis_r, nees_r, nmee_r, confidence] = confidence_interval(Nx, Nz, Nmonte)

alpha   = 0.05;

nees_r1 = chi2inv(alpha/2, Nmonte*Nx) /Nmonte;
nees_r2 = chi2inv(1-(alpha/2), Nmonte*Nx) /Nmonte;
nees_r  = [nees_r1; nees_r2]; 

nis_r1  = chi2inv(alpha/2, Nmonte*Nz)/Nmonte;
nis_r2  = chi2inv(1-(alpha/2), Nmonte*Nz)/Nmonte;

nis_r   = [nis_r1; nis_r2]; 

nmee_r  = 1.96/sqrt(Nmonte);
nmee_r  = [nmee_r; -nmee_r]; 

confidence = 100*(1-alpha);

end 
% .........................................................................






