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
