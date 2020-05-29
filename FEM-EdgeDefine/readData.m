function [dn,U,F,mp,nodemap,nd,ndim,ndn,ne,nen,nl,nm,nn,np,Gcoord,TH,nq,LC,NGP] = readData(finp)
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
area = X*Y;
%--------------------
DUMMY = fgets(finp);
TMP = str2num(fgets(finp));
[nmpc,np,NGP] = deal(TMP(1),TMP(2),TMP(3));
DUMMY = fgets(finp);
TMP = str2num(fgets(finp));
[seedx,seedy,LM] = deal(TMP(1),TMP(2),TMP(3));
ne = (seedx-1)*(seedy-1);
% ============================================================
% Global coordinate
Gcoord = zeros(seedx*seedy,2);
for ix = 1:seedx-1
    x_coord(ix+1) = ix*X/(seedx-1);
end
for iy = 1:seedy-1
    y_coord(iy+1) = iy*Y/(seedy-1);
end
nn = seedx * seedy; inn = 1:1:nn;
for j = 1:seedy
    for i  = 1:seedx
        
        n = i+seedx*(j-1);
        jseedx = x_coord(i);iseedy = y_coord(j);
        Gcoord(inn(n),:) = [jseedx, iseedy];
    end
end
% for in = 1:nn
%     Gcoord(in,1) = in;
% end
%============================================================
%Global Mesh Mapping
[nn,ne,nodemap,Gcoord] = MeshGenerate(seedx,seedy,nen,ne,nn,Gcoord,LM);
%============================================================
% =====================================================
% Local Mesh Mapping

if LM == 0
    disp('No local mesh refinement');
else
    DUMMY = fgets(finp);
    for i = 1:LM
        TMP = str2num(fgets(finp));
        %         [localrefineno(i,:),localx(i,:),localy(i,:)] = deal(TMP(1),TMP(2),TMP(3));
        [lfn,localx,localy] = deal(TMP(1),TMP(2),TMP(3));
        if i == 1
            lfn = lfn;
        else
            lfn = lfn-i+1;
        end
        lne = (localx-1)*(localy-1);lnn = localx*localy;
        [Gcoord,nodemap] = Localrefinemesh1(lfn,lne,nodemap,Gcoord,localx,localy,nen,ndim,nn);
        nn = nn - nen + lnn;
        ne = ne - 1 + lne;
    end
end
% =====================================================
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
nq = nn*ndn;
%-----------------------------
DUMMY = fgets(finp);DUMMY = fgets(finp);
U = zeros(nq,1);
TMP = str2num(fgets(finp));
[edge,Xno,Yno,value] = deal(TMP(1),TMP(2),TMP(3),TMP(4));
[U,dn,node] = edgedefine(edge,Xno,Yno,value,U,seedx,seedy);
[n,~]=size(node);
nd = n*2;
% 
% for i = 1:nd
%     TMP = str2num(fgets(finp));
%     [dn(i,:),U(dn(i),:)] = deal(TMP(1),TMP(2));
% end
DUMMY = fgets(finp);DUMMY = fgets(finp);
F = zeros(nq,1);
TMP = str2num(fgets(finp));
[edge,Xno,Yno,value] = deal(TMP(1),TMP(2),TMP(3),TMP(4));
[F,~,node] = edgedefine(edge,Xno,Yno,value,F,seedx,seedy);
[n,~]=size(node);
nl = n;
% 
% for i = 1:nl
%     TMP = str2num(fgets(finp));
%     [il,F(il,:)] = deal(TMP(1),TMP(2));
% end
%
%------------------
DUMMY = fgets(finp);
for i = 1:nm
    TMP = str2num(fgets(finp));
    [i,mp(i,:),TH] = deal(TMP(1),TMP(2:np),TMP(np+1));
end

% for i = 1:localno
%     TMP = str2num(fgets(finp));
%     [localrefineno(i,:),localx(i,:),localy(i,:)] = deal(TMP(1),TMP(2),TMP(3));
% end

% MAT = zeros(ne,1);
% MAT(:,1) = 1;
% TH = zeros(ne,1);
% TH(:,1) = 10;


fclose(finp);
% PlotMesh(Gcoord,nodemap,ne,nn,nen)
end