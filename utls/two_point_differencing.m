function [xk, zk, x0_hat, P0] = two_point_differencing(model, xk, zk, Ts, R)

% Two point differencing filter initialization for second-order
% models [Section 5.5.3]
%
% x0_hat - intial state vector
% P0     - intial covariance matrix
%
% model - {'CWNA', 'DWNA', 'CWPA', 'DWPA'}
% zk - first two measurement samples at k = 1 and 2
% Ts - sampling period in secs
% R  - noise covariance matrix


switch upper(model)
    case {'CWNA', 'DWNA'}
        [x0_hat, P0] = TPD(zk(:,1:2), Ts, R);

    case {'CWPA', 'DWPA'}
        [x0_hat, P0] = TPD(zk(:,1:2), Ts, R);
        x0_hat = [x0_hat; 0];
        P0 = blkdiag(P0, 0);
end

% remove the two-points from the true sytsem and measurement vectors 
xk     = xk(:, 3:length(xk)); 
zk     = zk(:, 3:length(zk)); 


end

function [x0_hat, P0] =  TPD(zk, Ts, R)
% two point differencing using first two measurements 
x0_hat = [zk(1); (zk(2)-zk(1))/Ts]; 
P0     = [1 -1/Ts; 
          -1/Ts 2/Ts^2].*(R)^2; 
end


