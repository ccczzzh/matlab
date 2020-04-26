function[K,k]=stiffness(nn,mm,ne,nc,MAT,x,mp,TH)
k = zeros(nn,mm);

for i = 1:ne
    n1 = nc(i,1);               % one node
    n2 = nc(i,2);               % another node
    n3 = MAT(i);                % Material group indicator
    L = abs(x(n2)-x(n1));       % ELEMENT displacement
    EAL = mp(n3,1)*AREA(i)/L;   % k(ij) = EA/L
    
    k(n1,1) = k(n1,1) + EAL;
    k(n2,1) = k(n2,2) + EAL;
    
    IR = n1;
    if IR > n2
        IR = n2;
    end
    IC = floor(abs(n2-n1)) + 1;
    k(IR,IC) = k(IR,IC)-EAL;
end

K = zeros(nn,nn);

for i = 1:nn
    j = i;

    K(i,j) = k(i,1);
    if j <4
    j= j+1;
    K(i,j) = k(i,2);
    end
    
end
    
    
    
    
    