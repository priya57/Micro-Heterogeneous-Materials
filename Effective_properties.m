global r mat


% Young's modulus for the matrix
E1 = 100;

% ratio of the inclusion Young's modulus over the matrix Young's modulus
E2 = 1000;

% Poisson's ratio of the matrix
ny1 = .3;
ny2 = .2;

% kapa for the matrix
k1 = E1 / ( 3 * ( 1 - 2 * ny1))%55.556; 

% kapa for the inclusino
k2 = E2 / ( 3 * ( 1 - 2 * ny2))%83333.33;

% miu for the matrix
miu1 = E1 /  ( 2 * ( 1 + ny1)) %41.667;

% miu for the inclusion
miu2 = E2 /  ( 2 * ( 1 + ny2))%38461.538;

% radius of the inclusion
r_in = 1;

% length and width of the matrix
x = 12; 
y = 12;

% volume fraction of the inclusion
c    = pi * r_in^2/x/y;


% Effective material properties
%[ C_eff, k_eff, miu_eff, alpha, beta] = get_eff( k1, k2, miu1, miu2, .022);

%E_eff = 9 * k_eff * miu_eff / ( 3 * k_eff + miu_eff);

%ny_eff = ( 3 * k_eff - 2 * miu_eff) / ( 6 * k_eff + 2 * miu_eff);

if 0
    disp('');
    disp('---------Effective material parameter---------');
    disp( sprintf('Poisson''s ratio  : ny_eff    = %8.4f ', ny_eff));
    disp( sprintf('Young''s modulus  : E_eff     = %8.4f ', E_eff));
    disp( sprintf('Bulk modulus     : kappa_eff = %8.4f ', k_eff));
    disp( sprintf('Lame constant    : my_eff    = %8.4f ', miu_eff));
    disp('----------------------------------------------');
end

% Apply Linear displacement BC, u = eps0 * x
if 0
eps_0 = [ .1, 0;
          0, .1];
uniform_strain( eps_0);
end

% restore data
if 0
    mat = mat_backup;
    is_homo = 0;
end

% run program

% take reaction forces, w/o homogenization
if 0
    F_org = r;
end

% substitute effetive material properties into the model
if 1
    
    mat_backup = mat;
    for i = 1:length(mat)
        mat(i).E = E_eff * ones( size( mat(i).E));
        mat(i).ny = ny_eff * ones( size( mat(i).ny));
    end
end

% run the program

% take reaction forces, with homogenization
if 0
    F_org_homo = r;
end

% plot reaction forces
if 0
    figure
    plot( F_org, '*'); 
    hold on
    plot( F_org_homo, 'linestyle', 'none', 'marker', 'o', 'color', 'r');
    hold off
    title ('reaction forces');
    legend( ' Without homogenization', 'with homogenization');
end