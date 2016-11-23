%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Group Project-5
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% material properties
c_i = 0.25;
c_m = 1-c_i;
E_m = 100;
ny_m = 0.3;
E_i = 100000;
ny_i = 0.2;
disp(c_m)
%Ellipsoidal inclusion geometry details

a1 = 10;    %major axis
b = 0.4;    %scale factor
a2 = b*a1;  %minior axis

% kapa for the matrix
Kappa_m = E_m / ( 3 * ( 1 - 2 * ny_m));

% kapa for the inclusino
Kappa_i = E_i / ( 3 * ( 1 - 2 * ny_i));

% miu for the matrix
mu_m = E_m /  ( 2 * ( 1 + ny_m));

disp(mu_m);
% miu for the inclusion
mu_i = E_i /  ( 2 * ( 1 + ny_i));
disp(mu_i);

% calculation of Eshelby Tensor
S_i = Eshelby(a1,a2,mu_i)
S_m = Eshelby(a1,a2,mu_m)

C_i = mat_prop( Kappa_i, mu_i);
C_m = mat_prop( Kappa_m, mu_m);
disp('C_m')
disp(C_m)
disp('C_i')
disp(C_i)
h_C_Eff_1 = figure;
h_C_Eff_2 = figure;

% Effective Elasticity Tensor Voigt
C_Eff_V = zeros(6,6);
C_Eff_R = zeros(6,6);

if 1
    c_v_1 = zeros();
    c_v_2 = zeros();
    for l = 1:length(c_i)
        C_Eff_V = c_i(l)*C_i + (1-c_i(l))*C_m;
        disp('Eff_Voigt');
        disp(C_Eff_V);
        c_v_1(l) = C_Eff_V(1,1);
        c_v_2(l) = C_Eff_V(6,6);

    end
    % plot
    set( 0, 'currentfigure', h_C_Eff_1);
    hold on;
    h_C_Eff_V_1 = plot( c_i, c_v_1, 'g');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution');
    hold off;
    
    set( 0, 'currentfigure', h_C_Eff_2);
    hold on;
    h_C_Eff_V_2 = plot( c_i, c_v_2, 'g');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution');
    hold off;
end

% Effective Elasticity Tensor Reuss
if 0
    c_r_1 = zeros();
    c_r_2 = zeros();
    for l = 1:length(c_i)
        C_Eff_R = (c_i(l)*C_i^-1 + (1-c_i(l))*C_m^-1)^-1;
        disp('Eff_Reuss');
        disp(C_Eff_R);
        % plot

        c_r_1(l) = C_Eff_R(1,1);
        c_r_2(l) = C_Eff_R(6,6);

    end
    set( 0, 'currentfigure', h_C_Eff_1);
    hold on;
    h_C_Eff_R_1 = plot( c_i, c_r_1, '-r');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution');
    hold off;
    set( 0, 'currentfigure', h_C_Eff_2);
    hold on;
    h_C_Eff_R_2 = plot( c_i, c_r_2, '-r');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution');
    hold off;
end

% Effective Elasticity Tensor Dilute Distribution

if 0
    C_Eff_DD = zeros(6,6);
    c_dd_1 = zeros();
    c_dd_2 = zeros();
    for i = 1:length(c_i)
        [ C_Eff_DD] = dilute_analy_elli( c_i(i), Kappa_m, mu_m, Kappa_i, mu_i, ...
            'plane strain');
        disp('Eff_DD')
        disp(C_Eff_DD)
        c_dd_1(i) = C_Eff_DD(1,1);
        c_dd_2(i) = C_Eff_DD(6,6);
    end
    % plot
    set( 0, 'currentfigure', h_C_Eff_1);
    hold on;
    h1_kappa_DD_1 = plot( c_i, c_dd_1, '-b');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution', 'Location', 'NorthWest');
    hold off; 
    
    set( 0, 'currentfigure', h_C_Eff_2);
    hold on;
    h1_kappa_DD_2 = plot( c_i, c_dd_2, '-b');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution', 'Location', 'NorthWest');
    hold off;
end
    
if 0
    kappa_DF = zeros( length( c_i), 1);
    mu_DF    = zeros( length( c_i), 1);
    c_ds_1 = zeros();
    c_ds_2 = zeros();

    % c_i = 0 and c_i = 1 can't be computed by the Matlab server
    %kappa_DF(1) = kappa_m; kappa_DF(end) = kappa_i;
    %mu_DF(1)    = mu_m;    mu_DF(end)    = mu_i;

    for i = 1:length(c_i)
        [ kappa_DF(i), mu_DF(i)] = diff_analy( c_i(i), Kappa_m, mu_m, Kappa_i, mu_i, ...
            'iso');
        C_Eff_DS = mat_prop( kappa_DF(i), mu_DF(i));
        disp('C_Eff_DS');
        disp(C_Eff_DS);
        c_ds_1(i) = C_Eff_DS(1,1);
        c_ds_2(i) = C_Eff_DS(6,6);
    end
    % plot
    set( 0, 'currentfigure', h_C_Eff_1);
    hold on;
    h1_kappa_DS_1 = plot( c_i, c_ds_1, '-y');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution','differential Scheme', 'Location', 'NorthWest');
    hold off;
    
    set( 0, 'currentfigure', h_C_Eff_2);
    hold on;
    h1_kappa_DS_2 = plot( c_i, c_ds_2, '-y');
    legend( 'Voigt', 'Reuss', 'Dilute Distribution','differential Scheme', 'Location', 'NorthWest');
    hold off;
end
