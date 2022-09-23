function [nis_bounds, nees_bounds, nmee_bounds, confidence] = confidence_bounds(Nx, Nz, Nmonte)

alpha = 0.05;
confidence = 100*(1-alpha);

nees_r1 = chi2inv(alpha/2, Nmonte*Nx) /Nmonte;
nees_r2 = chi2inv(1-(alpha/2), Nmonte*Nx) /Nmonte;
nees_bounds = [nees_r1; nees_r2]; 

nis_r1 = chi2inv(alpha/2, Nmonte*Nz)/Nmonte;
nis_r2 = chi2inv(1-(alpha/2), Nmonte*Nz)/Nmonte;

nis_bounds = [nis_r1; nis_r2]; 

nmee_bounds = 1.96/sqrt(Nmonte);
nmee_bounds = [nmee_bounds; -nmee_bounds]; 

end 