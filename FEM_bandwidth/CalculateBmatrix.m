function [B,Jacob] = CalculateBmatrix(nen,node,ndim,nodecoord,NGP,GQ,ne)

for ine = 1:ne
    B = zeros(2,nen*ndim); % for 3D problem b = zeros(6,nq)
    local_coord = nodecoord(node(ine,:),:);
    for iNGP = 1:NGP
        
        [gdN,Jacob] = calcShapefunctionAndJacob(iNGP,nen,GQ.point,local_coord);
        
        factor = det(Jacob)*GQ.weight(iNGP);
        
        for iinen = 1:nen
            i1 = (iinen-1)*ndim +1;
            i2 = i1 +1;
            B(1,i1) = gdN(1,iinen);
            B(2,i2) = gdN(2,iinen);
            B(3,i1) = gdN(2,iinen);
            B(3,i2) = gdN(1,iinen);
        end
    end
end