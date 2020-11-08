function [D] = calcDValue(ne,MAT,mp,LC,TH)
for i = 1:ne
    M = MAT(i);
    E = mp(M,1);
    rho = mp(M,2);
    alpha = mp(M,3);
    %--- D() Matrix
    if LC == 1
        C1 = (E*TH(i))/(1-rho^2);
        C2 = C1*rho;
    elseif LC == 2
        C = (E*TH(i))/((1+rho)*(1-2*rho));
        C1 = C*(1-rho);
        C2 = C * rho;
    end
    C3 = (E*TH(i))*0.5/(1+rho);
    
    D(1,1) = C1;D(1,2) = C2;D(1,3) = 0;
    D(2,1) = C2;D(2,2) = C1;D(2,3) = 0;
    D(3,1) = 0;D(3,2) = 0;D(3,3) = C3;
end