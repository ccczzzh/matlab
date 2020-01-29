function [nfdf,idd,idn,idr,LM0,LM] = generateID(numn,nndf,numel,nen,...
                                     elc,vdbc,vnbc,vrbc)

% Initialize id array
idd  = zeros(numn,nndf);
idn  = zeros(numn,nndf);
idr  = zeros(numn,nndf);
nfdf = 0;

for ii = 1:numn
    for jj = 1:nndf
    N1 = vdbc(:,1) == ii;
    N2 = vnbc(:,1) == ii;
    N3 = vrbc(:,1) == ii;
    if vdbc(N1,jj+1) == 1
        idd(ii,jj) = -1;
    else
        nfdf = nfdf + 1;
        idd(ii,jj) = nfdf;
    end
    if vnbc(N2,jj+1) == 1
        idn(ii,jj) = -1;
    end
    if vrbc(N3,jj+1) == 1
        idr(ii,jj) = -1;
    end
    end
end

% generate LM0 and LM array
LM0 = zeros(nndf*nen,numel);
for a=1:nen
    for b = 1:nndf
        LM0(2*(a-1)+b,:) = elc(a,:).*2 - (nndf-b);
    end
end

LM = zeros(nndf*nen,numel);
for ee = 1:numel
    for aa = 1:nen
        for ii = 1:nndf
            p = (aa - 1)*nndf + ii;
            LM(p,ee) = idd(elc(aa,ee),ii);
        end
    end
end

end   % end of function generateID