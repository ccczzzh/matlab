function[gdN,Jacob] = calcShapefunctionAndJacob(k,nen,GQpoint,e_coord)
if nen == 3
    ksi = GQpoint(k,1);
    eta = GQpoint(k,2);
    N(1) = 1.0 - ksi - eta;
    N(2) = ksi;
    N(3) = eta;
    
    % ksi derivatives
    dN(1,1) = -1.0;dN(1,2) = 1.0;dN(1,3) = 0.0;
    % eta derivatives
    dN(2,1) = -1.0;dN(2,2) = 0.0;dN(2,3)= 1.0;
elseif nen == 4
    ksi = GQpoint(k,1);
    eta = GQpoint(k,2);
    
    % ksi derivatives
    dN(1,1) = -0.25*(1-eta);
    dN(1,2) = -dN(1,1);
    dN(1,3) = 0.25*(1+eta);
    dN(1,4) = -dN(1,3);
    
    % eta derivatives
    dN(2,1) = -0.25*(1-ksi);
    dN(2,2) = -0.25*(1+ksi);
    dN(2,3) = -dN(2,2);
    dN(2,4) = -dN(2,1);
    
    N(1) = 0.25*(1+eta-ksi+ksi*eta);
    N(2) = 0.25*(1+eta-ksi-ksi*eta);
    N(3) = 0.25*(1-eta-ksi+ksi*eta);
    N(4) = 0.25*(1-eta+ksi-ksi*eta);
end

%Calculate the Jacobian Matrix
Jacob = dN*e_coord;
gdN = Jacob\ dN;
end