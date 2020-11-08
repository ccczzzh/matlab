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
finp = fopen(fid1,'r');
tic;
[dn,nn,U,F,mp,node,ndim,ne,nen,nodecoord,TH,MAT,nq,LC,NGP,ndn] = readData(finp);
% Original X Y coordinate
UX = nodecoord(:,1);UY = nodecoord(:,2);
if ndim == 2
    UT = sqrt(UX.^2+UY.^2);
elseif ndim == 3
    UZ = U(:,3);
    UT = sqrt(UX.^2+UY.^2+UZ.^2);
end
[D] = calcDValue(ne,MAT,mp,LC,TH);
GQ = calcGaussPoint(nen,NGP);
[K,B,Jacob] = Kmat(nq,nen,node,ndim,nodecoord,NGP,GQ,D,ne);


DB = zeros(3,6);
for i = 1:3
    for j = 1:6
        for k = 1:3
            DB(i,j) = DB(i,j) + D(i,k) * B(k,j);
        end
    end
end

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
c = 0;
for i = 1:ne
    I1 = node(i,1);
    I2 = node(i,2);
    I3 = node(i,3);
    Q(1) = U(I1*2-1);
    Q(2) = U(I1*2);
    Q(3) = U(I2*2-1);
    Q(4) = U(I2*2);
    Q(5) = U(I3*2-1);
    Q(6) = U(I3*2);
    for j = 1:3
        for k = 1:6
            c = c + DB(j,k) * Q(k);
        end
        STRESS(i,j) = c;
    end
    S1 = STRESS(i,1);
    S2 = STRESS(i,2);
    S3 = STRESS(i,3);
    if S3 == 0
        if S1 > S2
            Ang = 0;
        elseif S2 > S1
            Ang = 90;
        end
    else
        A = 0.5 * (S1 + S2);
        B = sqrt((0.5 * (S1 - S2))^2 + S3^2);
        Ang = 0.5 * atan((2*S3)/(S1-S2));
        PS1 = A + B;
        PS2 = A - B;
    end
    PRINSTRESS(i,1) = PS1;
    PRINSTRESS(i,2) = PS2;
    PRINSTRESS(i,3) = Ang;
    vms(i) = sqrt(STRESS(i,1)^2 - STRESS(i,1) * STRESS(i,2) - STRESS(i,2)^2 +3*STRESS(i,3)^2);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% count = 1;
% for i = 1:ne
%     
%     for j = 1:nen
%         if j == nen
%             x2 = (UX(node(i,j))-UX(node(i,1)))^2;
%             y2 = (UY(node(i,j))-UY(node(i,1)))^2;
%             x21 = (UXnew(node(i,j))-UXnew(node(i,1)))^2;
%             y21 = (UYnew(node(i,j))-UYnew(node(i,1)))^2;
%         else
%             x2 = (UX(node(i,j))-UX(node(i,(j+1))))^2;
%             y2 = (UY(node(i,j))-UY(node(i,(j+1))))^2;
%             x21 = (UXnew(node(i,j))-UXnew(node(i,(j+1))))^2;
%             y21 = (UYnew(node(i,j))-UYnew(node(i,(j+1))))^2;
%         end
%         original_length(count) = sqrt(x2+y2); 
%         delta_length(count) = sqrt(x21+y21);
%         count = count + 1;
%     end
% end
% epsilone = delta_length./original_length;
% S = mp(1).*epsilone;
% 
for i = 1:ne
   meshmaping(i,:) = [node(i,1) node(i,2) node(i,3) node(i,1)];
end

% 
% node = node';
% count = 1;
% for j = 1 : ne
%     for i =1:nen
%        inode(count) = node(i,j); 
%        count= count+1;
%     end
% end
% [~,n]=size(S);
% nodeStress = zeros(1,nn);
% count = 1;
% 
% for j = 1:nn
%     for i = 1:n
%         if inode(i) == j
%             nodeStress(j) = nodeStress(j) + S(i);
%             continue
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

node = node';

%-----------------------------------------------------
%                       Plot Mesh
%-----------------------------------------------------
PlotMesh(nodecoord,node,ne,nn,nen);
Unew = nodecoord + Unew;
PlotdeformMesh(Unew,node,ne,nen);

UXnew = Unew(:,1);UYnew = Unew(:,2);
for i = 1:ne
   meshX(i,:) = UXnew(meshmaping(i,:));
   meshY(i,:) = UYnew(meshmaping(i,:));
   meshStress(i,:) = vms(meshmaping(i,:));
end
figure;hold on;
for i = 1:ne
    fill(meshX(i,:),meshY(i,:),meshStress(i,:))
    
end

hold off
colorbar;
title('von Mises Stress');
axis equal
fprintf(1,'\nCalculated unknowns are \n\n');
fprintf(' Node              x                  y                   ux                    uy \n');
fprintf(1,'======================================================================================\n');

for i = 1:nn
    fprintf(1,' %-5d %18.8f %18.8f %20.8f %20.8f\n', i,UX(i),UY(i),UXnew(i),UYnew(i));
end
toc;
