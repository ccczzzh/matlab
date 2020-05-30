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
tic;
% fid1 = input('Input Data File Name: ','s');
% finp = fopen(fid1,'r');
finp = fopen('FEM2D_mesh.inp');
% fid2 = input('Input Data File Name: ','s');
% foup = fopen(fid2,'w');

[dn,U,F,mp,node,nd,ndim,ndn,ne,nen,nl,nm,nn,np,nodecoord,TH,nq,LC,NGP]...
    = readData(finp);
UX = nodecoord(:,1);UY = nodecoord(:,2);% UZ = U(:,3);
UT = sqrt(UX.^2+UY.^2);

%[mm] = MatrixModification(node,ne);
% Calculate the Gauss point
GQ = calcGaussPoint(nen,NGP);
[D] = calcDValue(mp,TH,LC);

[K,DB] = Kmat(nq,nen,node,ndim,nodecoord,NGP,GQ,D,ne);
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
for i = 1:ne
%     [DA, K,DB,kmatrix,B,D] = dbmat(ne, LC, mp, node, nodecoord,nen,nq);
    E = mp(1);rho = mp(2);alpha = 0;
    % 4-node elements
    I1 = node(i,1);I2 = node(i,2);I3 = node(i,3);I4 = node(i,4);
    Q(1) = F(2 * I1 - 1);
    Q(2) = F(2 * I1);
    Q(3) = F(2 * I2 - 1);
    Q(4) = F(2 * I2);
    Q(5) = F(2 * I3 - 1);
    Q(6) = F(2 * I3);
    Q(7) = F(2 * I4 - 1);
    Q(8) = F(2 * I4);
    C1 = alpha * TH;
    if LC == 2
        C1 = C1 * (1 + rho);
    end
    for ii = 1:3
        C = 0;
        for k = 1:8
            C = C + DB(ii,k) * Q(k);
        end
        STR(ii) = C - C1 * (D(ii,1)+D(ii,2));
    end
    if STR(3) == 0
        S1 = STR(1);S2 = STR(2);
        ANG = 0;
        if S2 > S1
            S1 = STR(2);S2 = STR(1);
            ANG = 90;
        end
    else
        C = 0.5 * (STR(1) + STR(2));
        R = sqrt(0.25 * (STR(1) - STR(2))^2 + STR(3) ^2);
        S1 = C + R;
        S2 = C - R;
        if C > STR(1)
            ANG = 57.2957795 * atan(STR(3) / (STR(3))^2);
            if STR(3) >0;ANG = 90 - ANG; end
            if STR(3) >0;ANG = -90-ANG;end
        else
            ANG = 57.29577951 * atan(STR(3) / (STR(3) / (STR(1)-S2)));
        end
    end
    Stress(i,1) = STR(1);
    Stress(i,2) = STR(2);
    Stress(i,3) = STR(3);
    PrincipalStress(i,1) = S1;
    PrincipalStress(i,2) = S2;
    PrincipalStress(i,3) = ANG;
end

%-----------------------------------------------------
%                       Plot Mesh
%-----------------------------------------------------
PlotMesh(nodecoord,node,ne,nn,nen);
D_U = Unew;
Unew = nodecoord + Unew;
PlotdeformMesh(Unew,node,ne,nen);
strain = D_U./Unew;
stress = mp(1).*strain;
%----------------------------------------------------
%       Table Print section
%----------------------------------------------------

fprintf(1,'\nCalculated unknowns are \n\n');
fprintf(' Node              x                  y                   ux                    uy \n');
fprintf(1,'======================================================================================\n');

for i = 1:nn
    fprintf(1,' %-5d %18.8f %18.8f %20.8f %20.8f\n', i,UX(i),UY(i),UXnew(i),UYnew(i));
end
toc;
