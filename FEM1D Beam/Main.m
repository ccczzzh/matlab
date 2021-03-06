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
%global nn ne nm ndim nen ndn
fid1 = input('Input Data File Name: ','s');
finp = fopen(fid1,'r');
% fid2 = input('Input Data File Name: ','s');
% foup = fopen(fid2,'w');
%====================================================
% INPUT STEP
[dn,U,F,mp,nc,nd,ndim,ndn,ne,nen,nl,nm,nn,np,x,AREA,MAT]...
    = readData(finp);
[mm] = MatrixModification(nc,ne);

% Stiffness Matrix
[K,k]=stiffness(nn,mm,ne,nc,MAT,x,mp,AREA);

%----------------------------------------------
%       Calculate Stiffness Matrix
%----------------------------------------------
move = setdiff(1:nn,dn);
Kreduced = K(move,move);
Freduced = F(move);
Ureduced = Kreduced\Freduced;
U(move) = Ureduced;
F = K*U;
%------------------------------------------
%        Stress Calculation
%------------------------------------------
for n = 1:ne
    n1 = nc(n,1);n2 = nc(n,2);n3=MAT(n);
    PR = (F(n2)-F(n1))/(U(n2)-U(n1));
    S(n) = mp(n3,1)*PR;  
end

