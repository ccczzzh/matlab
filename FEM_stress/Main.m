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
fid1 = input('Input Data File Name: ','s');
% finp = fopen('FEM2D.inp');
finp = fopen(fid1,'r');
% fid2 = input('Input Data File Name: ','s');
% foup = fopen(fid2,'w');
tic;
[dn,nn,U,F,mp,node,ndim,ne,nen,nodecoord,TH,MAT,nq,LC,NGP] = readData(finp);
% Original X Y coordinate
UX = nodecoord(:,1);UY = nodecoord(:,2);% UZ = U(:,3);
UT = sqrt(UX.^2+UY.^2);
%UT = sqrt(UX.^2+UY.^2+UZ.^2);


GQ = calcGaussPoint(nen,NGP);
[D] = calcDValue(ne,MAT,mp,LC,TH);

[K,B] = Kmat(nq,nen,node,ndim,nodecoord,NGP,GQ,D,ne);


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
% Deformed X Y Coordinate
UXnew = Unew(:,1);UYnew = Unew(:,2);
% %------------------------------------------
% %        Stress Calculation
% %------------------------------------------
% % % % for n = 1:ne
% % % %     n1 = nc(n,1);n2 = nc(n,2);n3=MAT(n);
% % % %     PR = (F(n2)-F(n1))/(U(n2)-U(n1));
% % % %     S(n) = mp(n3,1)*PR;  
% % % % end
% % % % count = 1;
% % % % for i = 1:ne
% % % %     
% % % %     for j = 1:nen
% % % %         if j == nen
% % % %             x2 = (UX(node(i,j))-UX(node(i,1)))^2;
% % % %             y2 = (UY(node(i,j))-UY(node(i,1)))^2;
% % % %             x21 = (UXnew(node(i,j))-UXnew(node(i,1)))^2;
% % % %             y21 = (UYnew(node(i,j))-UYnew(node(i,1)))^2;
% % % %         else
% % % %             x2 = (UX(node(i,j))-UX(node(i,(j+1))))^2;
% % % %             y2 = (UY(node(i,j))-UY(node(i,(j+1))))^2;
% % % %             x21 = (UXnew(node(i,j))-UXnew(node(i,(j+1))))^2;
% % % %             y21 = (UYnew(node(i,j))-UYnew(node(i,(j+1))))^2;
% % % %         end
% % % %         original_length(count) = sqrt(x2+y2);
% % % %         delta_length(count) = sqrt(x21+y21);
% % % %         count = count + 1;
% % % %     end
% % % % end
% % % % epsilone = delta_length./original_length;
% % % % S = mp(1).*epsilone;

% e_strain = D * B %.* U
count = 1;
local_U = zeros(nen*2,ne);
for i = 1:ne
    local_ux(:,i) = UXnew(node(i,:));
    local_uy(:,i) = UYnew(node(i,:));
    local_U(1:3,i) = local_ux(:,i);
    local_U(4:6,i) = local_uy(:,i);

%     local_ux(1:count*3,i) = U(location1);
%     
%     local_uy(1:count*3,i) = U(location2);
%     count = count +1;
end
Stress = D * B * local_U
%-----------------------------------------------------
%                       Plot Mesh
%-----------------------------------------------------
PlotMesh(nodecoord,node,ne,nn,nen);
Unew = nodecoord + Unew;
PlotdeformMesh(Unew,node,ne,nen);
fprintf(1,'\nCalculated unknowns are \n\n');
fprintf(' Node              x                  y                   ux                    uy \n');
fprintf(1,'======================================================================================\n');

for i = 1:nn
    fprintf(1,' %-5d %18.8f %18.8f %20.8f %20.8f\n', i,UX(i),UY(i),UXnew(i),UYnew(i));
end
toc;
