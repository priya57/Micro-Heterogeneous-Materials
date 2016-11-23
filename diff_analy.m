function [ kappa_eff, mu_eff] = diff_analy(c_i, Kappa_m, mu_m, Kappa_i, mu_i, str)
% c_i           = volume fraction of inclusion
% kappa_m, mu_m = kappa, mu of matrix phase
% kappa_i, mu_i = kappa, mu of inclusion phase
switch str
    case 'iso'
        z0       = [Kappa_m,mu_m] ; % Initial conditions for the ODE
        v2_range = [0,c_i] ;          % Range of integration

        % Standard MATLAB ODE integrator
        [v2_i,zi] = ode45(@DFode,v2_range,z0,[],Kappa_i,mu_i);

        % Values at the end of the integration
        z(1:2) = zi(length(zi),1:2) ;

        kappa_eff = z(1);
        mu_eff    = z(2);
    case 'plain strain'
end


% DF MODEL Ordinary Differential Equation
function dfdv2 = DFode(c_i,eff,Kappa_i,mu_i)
effk = eff(1) ;
effu = eff(2) ;
% right hand side of ODE (output must be a column vector)
dfdv2(1,1) = 1/(1-c_i)*(Kappa_i-effk)*((3*effk+4*effu)/(3*Kappa_i+4*effu));
dfdv2(2,1) = 1/(1-c_i)*(mu_i-effu)*((5*effu*(3*effk+4*effu))/(effu*(9*effk+8*effu)+6*mu_i*(effk+2*effu)));

