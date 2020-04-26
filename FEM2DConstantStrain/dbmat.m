function [DA, K,DB] = dbmat(ne, LC, MAT, mp, nc, x,K,kc,cn)
%----- D(), B() and DB() matrices

%  --- Material Properties
    for i = 1:ne
        M = MAT(i);
        E = mp(M,1);
        rho = mp(M,2);
        alpha = mp(M,3);
%--- D() Matrix 
    if(LC == 1)
        C1 = E/(1-rho^2);
        C2 = C1*rho;
    else
        C = E/((1+rho)*(1-2*rho));
        C1 = C*(1-rho);
        C2 = C * rho;
    end
     C3 = E*0.5/(1+rho);
     
     D(1,1) = C1;
     D(1,2) = C2;
     D(1,3) = 0;
     D(2,1) = C2;
     D(2,2) = C1;
     D(2,3) = 0;
     D(3,1) = 0;
     D(3,2) = 0;
     D(3,3) = C3;
 %--- Strain-Displacement Matrix B()
     I1 = nc(i,1);
     I2 = nc(i,2);
     I3 = nc(i,3);
     X1 = x(I1,1);
     Y1 = x(I1,2);
     X2 = x(I2,1);
     Y2 = x(I2,2);
     X3 = x(I3,1);
     Y3 = x(I3,2);
     X32 = X3-X2;
     X13 = X1-X3;
     X21 = X2-X1;
     Y23 = Y2-Y3;
     Y31 = Y3-Y1;
     Y12 = Y1-Y2;
     DA = X13*Y23-X32*Y31;
     % Definition of B matrix
     B(1,1) = Y23/DA;
     B(2,1) = 0;
     B(3,1) = X32/DA;
     B(2,1) = 0;
     B(2,2) = X32/DA;
     B(3,2) = Y23/DA;
     B(1,3) = Y31/DA;
     B(2,3) = 0;
     B(3,3) = X13/DA;
     B(1,4) = 0;
     B(2,4) = X13/DA;
     B(3,4) = Y31/DA;
     B(1,5) = Y12/DA;
     B(2,5) = 0;
     B(3,5) = X21/DA;
     B(1,6) = 0;
     B(2,6) = X21/DA;
     B(3,6) = Y12/DA;
     % DB matrix
     for a = 1:3
         for b = 1:6
             C = 0;
             for c = 1:3
                 C = C + D(a, c) * B(c, b);
             end
             DB(a, b) = C;
         end
     end
     
     for a = 1:6
         for b = 1:6
             C = 0;
             for c = 1:3
                 C = C + B(c,a)*DB(c,b);
             end
             kmatrix(a,b) = C;
         end
     end

     for a = 1:6
         kc(2*a-1) =cn(a)*2-1;
         kc(2*a) = cn(a)*2;
     end

     for a = 1:6
         for b = 1:6
             K(kc(a,i),kc(b,i)) = K(kc(a,i),kc(b,i))+kmatrix(a,b);
         end
     end
    end
end
