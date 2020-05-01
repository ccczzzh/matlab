function [dn,U,F,mp,nodemap,nd,ndim,ndn,ne,nen,nl,nm,nn,np,coord,TH,MAT,nq,LC,NGP] = readData(finp)
% disp('  1) Plane Stress Analysis');
% disp('  2) Plane Strain Analysis');
% LC = input('  Choose 1(default) or 2 :');
% if isempty(LC)| LC < 1 | LC > 2
%     LC = 1;
% end
LC = 1;
TITLE = fgets(finp);
DUMMY = fgets(finp);
TMP = str2num(fgets(finp));
[X,Y,nm,ndim,nen,ndn] = deal(TMP(1),TMP(2),TMP(3),...
    TMP(4),TMP(5),TMP(6));

%--------------------
DUMMY = fgets(finp);
TMP = str2num(fgets(finp));
[nd,nl,nmpc,np,NGP] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5));
DUMMY = fgets(finp);
TMP = str2num(fgets(finp));
[SeedX,SeedY] = deal(TMP(1),TMP(2));
%-------------------
% DUMMY = fgets(finp);
% for i = 1:nn
%     TMP = str2num(fgets(finp));
%     [i,coord(i,:)] = deal(TMP(1),TMP(2:1+ndim));
% end
% %------------------
% DUMMY = fgets(finp);
% for i = 1:ne
%     TMP = str2num(fgets(finp));
%     [i,nodemap(i,:),MAT(i,:),TH(i,:)] = ...
%         deal(TMP(1),TMP(2:1+nen),TMP(2+nen),TMP(3+nen));
% end
%------------------
[nn,ne,nodemap,coord] = MeshGenerate(X,Y,SeedX,SeedY,nen);
MAT = zeros(ne,1);
MAT(:,1) = 1;
TH = zeros(ne,1);
TH(:,1) = 10;
nq = nn*ndn;
%-----------------------------
DUMMY = fgets(finp);
U = zeros(nq,1);
for i = 1:nd
    TMP = str2num(fgets(finp));
    [dn(i,:),U(i,:)] = deal(TMP(1),TMP(2));
end
DUMMY = fgets(finp);
F = zeros(nq,1);
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