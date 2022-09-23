function [xk_hat, Pk, nis] = kalman_filter(xkm1_hat, Pkm1, zk, F, Q, H, R)
% Kalman filter [Section 5.2.4]

% xk_hat - current estimate at time index k
% Pk     - current covariance estimate
% NIS    - normalized innovation error squared
%
% xkm1_hat - previous estimate at time index k-1
% Pkm1     - previous covariance estimate at time index k-1
% zk       - current measurement at time index k

% allocate memory
Nsamples = size(zk,2);
Nx       = size(xkm1_hat,1);
xk_hat   = zeros(Nx, Nsamples);
Pk       = zeros(Nx, Nx, Nsamples);
nis      = zeros(1, Nsamples);

if Nsamples == 1 % one run filter
    [xk_hat, Pk, nis] = KF(xkm1_hat, Pkm1, zk(:,1), F, Q, H, R);
else % batch filter
    for k = 1:Nsamples
        if k == 1 % intialize the filter
            xk_hat(:,1) = xkm1_hat;
            Pk(:,:,1)   = Pkm1;
        else
            [xk_hat(:,k), Pk(:,:,k), nis(:,k)] = KF(xk_hat(:,k-1), Pk(:,:,k-1), zk(:,k), F, Q, H, R);
        end
    end
end
end

function [xk_hat, Pk, nis] = KF(xkm1_hat, Pkm1, zk, F, Q, H, R)
% State and covariance prediction
xpred  = F*xkm1_hat;    % predicted state
Ppred  = F*Pkm1*F' + Q; % predicted covariance
% Measurement prediction
zk_hat = H*xpred;             % predicted measurement
inov   = zk - zk_hat;         % residual/innovation
S      = R + H*Ppred*H';      % innovation covariance
nis    = inov'*S^(-1)*inov;   % normalized innovation squared
% State and covariance update
G      = Ppred*H'/(S);     % kalman gain
xk_hat = xpred + G*inov;   % updated state
% Pk   = Ppred - G*S*G';   % updated covariance
Pk     = (eye(size(Q,1)) - G*H)*Ppred*(eye(size(Q,1)) - G*H)' + G*R*G'; % Joseph form (alternative covarian ce update for numerical stability)
end