function [ C_Eff_DD] = dilute_analy_elli( c_i, Kappa_m, mu_m, Kappa_i, mu_i, str)
%obtaining the analytical effective material properties by dilute
%distribution method
%ci:   volume fraction of the inclusion
%kappa_m:  kappa of the matrix material
%mu_m:     mu of the matrix material
%kappa_i:  kappa of the inclusion
%mu_i:     mu of the inclusion

switch str
    case 'iso'
        
    case 'plane strain'
        a1 = 10;    %major axis
        b = 0.4;    %scale factor
        a2 = b*a1;  %minior axis
        % kapa for the matrix
                
        C_i = mat_prop( Kappa_i, mu_i);
        C_m = mat_prop( Kappa_m, mu_m);
        S_m = Eshelby(a1,a2,mu_m);
        
        L = (eye(6,6)+S_m.*(1-C_i).^-1.*(C_i-C_m));
        C_Eff_DD = C_m + c_i * ( C_i - C_m) .* L;
        
end

