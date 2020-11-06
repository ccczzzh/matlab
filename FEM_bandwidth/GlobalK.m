function[Ks] =GlobalK(ndim,nn,BW,D,B,Jacob,node,ne,nen,ndn)
%============================================
nq = ndim*nn;
% Global Stiffness Matrix
Ks = zeros(nq,BW);
DB = zeros(3,6);
for i = 1:3
    for j = 1:6
        for k = 1:3
            DB(i,j) = DB(i,j) + D(i,k) * B(k,j);
        end
    end
end
for n = 1:ne
    for i = 1:6
        for j = 1:6
            Ktemp = 0;
            for k = 1:3
                Ktemp = Ktemp +0.5*abs(det(Jacob))*B(k,i)*DB(k,j);
            end
            KE(i,j)=Ktemp;
        end
    end
%     % Temperature Load Vector
%     M = MAT(N);
    % Global Location
    for ii = 1:nen
        NRT = ndn*(node(n,ii)-1);
        for it = 1:ndn
            NR = NRT +it;
            i = ndn *(ii-1) + it;
            for jj = 1:nen
                NCT = ndn*(node(n,jj)-1);
                for jt = 1:ndn
                    j = ndn *(jj-1)+ jt;
                    NC = NCT + jt-NR +1;
                    if NC>0
                        Ks(NR,NC) = Ks(NR,NC) +KE(i,j);
                    end
                end
            end
            % if temperature load included
           % F(NR) = F(NR) +TL(i);
        end
    end
end


%============================================
