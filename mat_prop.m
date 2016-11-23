function C = mat_prop(kappa, mu)
    C = zeros(6,6);
    C(1,1) = kappa+(4/3*mu);
    C(1,2) = kappa-(2/3*mu);
    C(4,4) = 2*mu;
    C(5,5) = C(4,4);
    C(2,2) = C(1,1);
    C(3,3) = C(2,2);
    C(6,6) = C(5,5);
    C(1,3) = C(1,2);
    C(2,1) = C(1,2);
    C(3,1) = C(1,3);
    C(2,3) = C(1,2);
    C(3,2) = C(2,3);
end







