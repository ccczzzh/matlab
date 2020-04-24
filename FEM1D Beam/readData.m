function [dn,U,F,mp,nc,nd,ndim,ndn,ne,nen,nl,nm,nn,np,x,AREA,MAT] = readData(finp)
%global nn ne nm ndim nen ndn
% fid1 = input('Input Data File Name: ','s');
% finp = fopen(fid1,'r');
% fid2 = input('Input Data File Name: ','s');
% foup = fopen(fid2,'w');
TITLE = fgets(finp);
DUMMY = fgets(finp);
TMP = str2num(fgets(finp));
[nn,ne,nm,ndim,nen,ndn] = deal(TMP(1),TMP(2),TMP(3),...
    TMP(4),TMP(5),TMP(6));
%--------------------
DUMMY = fgets(finp);
TMP = str2num(fgets(finp));
[nd,nl,nmpc,np] = deal(TMP(1),TMP(2),TMP(3),TMP(4));
%-------------------
DUMMY = fgets(finp);
for i = 1:nn
    TMP = str2num(fgets(finp));
    [i,x(i,:)] = deal(TMP(1),TMP(2:1+ndim));
end
%------------------
DUMMY = fgets(finp);
for i = 1:ne
    TMP = str2num(fgets(finp));
    [i,nc(i,:),MAT(i,:),AREA(i,:)] = ...
        deal(TMP(1),TMP(2:1+nen),TMP(2+nen),TMP(3+nen));
end
%------------------
DUMMY = fgets(finp);
U = zeros(nn,1);
for i = 1:nd
    TMP = str2num(fgets(finp));
    [dn(i,:),U(i,:)] = deal(TMP(1),TMP(2));
end
DUMMY = fgets(finp);
F = zeros(nn,1);
for i = 1:nl
    TMP = str2num(fgets(finp));
    [i,F(i,:)] = deal(TMP(1),TMP(2));
end

%------------------
DUMMY = fgets(finp);
for i = 1:nm
    TMP = str2num(fgets(finp));
    [i,mp(i,:)] = deal(TMP(1),TMP(2:1+np));
end
fclose(finp);

end