function [y] = BENR(yn,dt)

m=40;
A=0.5;
cD=1.;
ro=1.2;
g=9.81;

tolerance=1.e-8;
maxiter=100;
iter=0;
y=yn;

while(iter<maxiter)
    iter=iter+1;


%Residual
R= m*(y-yn)/dt+(0.5*ro*cD*A)*y^2-m*g;
%Derivative
dR= m/dt+y*(ro*cD*A);
  
if(abs(dR)<1.e-14) break; end;
if(abs(R)<tolerance) break; end;

%update NR guess
y=y-R/dR;

end;