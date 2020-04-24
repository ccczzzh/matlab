function [mm] = MatrixModification(nc,ne)
mm = 0;
for n = 1:ne
    nABS = floor(abs(nc(n,1)-nc(n,2)))+1;
    if (mm < nABS)
        mm = nABS;
    end
end

end