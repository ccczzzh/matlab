function [Ks,F] = solveBandwidth(nn,BW,Ks,F)
% Forward Elimination
for k = 1:nn-1
    nbk = nn-k+1;
    if (nn-k+1)>BW
        nbk = BW;
    end
    for i = k+1:nbk+k-1
        i1 = i-k+1;
        c = Ks(k,i1)/Ks(k,1);
        for j = i:nbk+k-1
            j1 = j-i+1;
            j2 = j-k+1;
            Ks(i,j1) = Ks(i,j1)-c*Ks(k,j2);
        end
        F(i) = F(i)-c*F(k);
    end
end
% Back Substitution
F(nn) = F(nn)/Ks(nn,1);
for ii = 1:nn-1
    i=nn-ii;
    nbi = nn-i+1;
    if (nn-i+1)>BW
        nbi = BW;
    end
    sum = 0;
    for j = 2:nbi
        sum = sum+Ks(i,j)*F(i+j-1);
    end
    F(i) = (F(i)-sum)/Ks(i,1);
end