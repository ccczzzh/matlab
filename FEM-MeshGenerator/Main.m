%%% ------------------------------------------------
% NN ---- Number of Nodes
% NE ---- Number of Elements
% NM ---- Number of Different Materials
% NDIM -- Number of Coordinates per Node(NDIM = 2 for two dimensional problems)
% NEN --- Number of Nodes per Element(NEN = 3 for CST Elem)
% NDN --- Number of Degrees of Freedom per Node (NDN = 2 for BEAM)
% ND ---- Number of Specified Displacement Degrees of Freedom
% NL ---- Number of Applied Component Loads
% NMPC -- Number of Multipoint Constraints
% N1,N2 - Node1,Node2. Connectivity table
% MAT# -- Material Flag
% NOC --- Node coordinate
% MP --- Material Property
% NP --- Number of Property is used
% mm --- matrix modification index
% TP --- Thermal Property
% TV --- Thermal Value
%%%-------------------------------------------------
clc;clear all;close all;
% fid1 = input('Input Data File Name: ','s');
% finp = fopen(fid1,'r');
finp = fopen('FEM2D_mesh.inp');
% fid2 = input('Input Data File Name: ','s');
% foup = fopen(fid2,'w');

[dn,U,F,mp,node,nd,ndim,ndn,ne,nen,nl,nm,nn,np,nodecoord,TH,MAT,nq,LC,NGP]...
    = readData(finp);
UX = nodecoord(:,1);UY = nodecoord(:,2);% UZ = U(:,3);
UT = sqrt(UX.^2+UY.^2);

%[mm] = MatrixModification(node,ne);
% Calculate the Gauss point
GQ = calcGaussPoint(nen,NGP);
[D] = calcDValue(ne,MAT,mp,LC);

[K] = Kmat(nq,nen,node,ndim,nodecoord,NGP,GQ,D,ne);
% end
% end
% 
% end
%[DA, K,DB] = dbmat(ne, LC, MAT, mp, node, nodecoord,nen,nq);

% %----------------------------------------------
% %       Calculate Stiffness Matrix
% %----------------------------------------------
move = setdiff(1:nq,dn);
Kreduced = K(move,move);
Freduced = F(move);
Ureduced = Kreduced\Freduced;
U(move) = Ureduced;
F = K*U;
Unew = reshape(U,[2,nn]);
Unew = Unew.';
UXnew = Unew(:,1);UYnew = Unew(:,2);
% %------------------------------------------
% %        Stress Calculation
% %------------------------------------------
% for n = 1:ne
%     n1 = nc(n,1);n2 = nc(n,2);n3=MAT(n);
%     PR = (F(n2)-F(n1))/(U(n2)-U(n1));
%     S(n) = mp(n3,1)*PR;  
% end

%-----------------------------------------------------
%                       Plot Mesh
%-----------------------------------------------------
PlotMesh(nodecoord,node,ne,nn,nen);
Unew = nodecoord + Unew;
PlotdeformMesh(Unew,node,ne,nen)

