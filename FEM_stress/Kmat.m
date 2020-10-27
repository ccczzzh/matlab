function [K,B] = Kmat(nq,nen,node,ndim,nodecoord,NGP,GQ,D,ne)

Kinitial = zeros(nq,nq);
K = Kinitial;
kc = zeros(nen*2,1);
cn = node.';
for ine = 1:ne
    ke = zeros(nen*ndim,nen*ndim); % local K matrix
    B = zeros(3,nen*ndim); % for 3D problem b = zeros(6,nq)
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
        ke = ke + B' * D * B * factor;
        
        for a = 1:nen
            kc(2*a-1) = cn(a,ine)*2-1;
            kc(2*a) = cn(a,ine)*2;
        end
        for a = 1:nen*2
            for b = 1:nen*2
                K(kc(a,1),kc(b,1)) = K(kc(a,1),kc(b,1))+ke(a,b);
            end
        end
    end
end
end






