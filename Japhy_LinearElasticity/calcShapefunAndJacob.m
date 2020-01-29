function [N,gdN,Jacob] = calcShapefunAndJacob(k,nen,elemType,GQpoint,e_coord)

% Calculates shape functions and their derivatives 
N(1,nen) = 0.0;
dN(2,nen) = 0.0;
if elemType == 1  % Triangular element
    if nen == 3
        ksi = GQpoint(k,1);
        eta = GQpoint(k,2);
        N(1) = 1.0 - ksi - eta;
        N(2) = ksi;
        N(3) = eta;
            
        % ksi derivatives
        dN(1,1) = -1.0;
        dN(1,2) =  1.0;
        dN(1,3) =  0.0;     
            
        % eta derivatives of S
        dN(2,1) = -1.0;
        dN(2,2) =  0.0;
        dN(2,3) =  1.0;
    else
        disp('ERROR: For triangular elements only 3-node is supported.');
    end
    
elseif elemType == 2  % Quadrilateral element
    if nen == 4
        ksi = GQpoint(k,1);
        eta = GQpoint(k,2);
     
        % ksi derivatives
        dN(1,1) = -0.25*(1-eta);
        dN(1,2) = -dN(1,1);
        dN(1,3) =  0.25*(1+eta);
        dN(1,4) = -dN(1,3);     
        
        % eta derivatives
        dN(2,1) = -0.25*(1-ksi);
        dN(2,2) = -0.25*(1+ksi);
        dN(2,3) = - dN(2,2);
        dN(2,4) =  -dN(2,1);
        
        N(1) = dN(1,2) - ksi*dN(1,3);
        N(2) = dN(1,2) + ksi*dN(1,2);
        N(3) = dN(1,3) + ksi*dN(1,3);
        N(4) = dN(1,3) - ksi*dN(1,2);
    else
        disp('ERROR: For quadrilateral elements only 4-node is supported.');
    end
end
    
% Calculte Jacobian matrix
Jacob = dN * e_coord;
gdN = Jacob\ dN;   