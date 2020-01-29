%**********************************************************************
% Program for Linear Elastic Finite Element Analysis of Solids under
% Plane Stress, Plane Strain and Axisymmetric Conditions with 3-noded
% Triangular Elements
%
% Dr. E. A. de Souza Neto
%
% Civil and Computational Engineering Centre
% Swansea University
%
% MODULE: EG-M23 Finite Element Computational Analysis
%
% HISTORY
% November 2003:   Initial coding
% January  2011:   Code improvements made. Most loops removed; Reaction
%                  vector plot added
%***********************************************************************
clear;hold off;
% Set some basic variables
% ========================
% Number of degrees of freedom per node
ndofn=2;
%
% Data input phase. Read information from data file
% =================================================
%
filename = input('Input data file name -------> ','s');
fid = fopen(filename, 'r');

title = fscanf(fid, 'TITLE = %s',1);
% 
nnode = fscanf(fid, '\nNUMBER_OF_NODES_PRE_ELEMENT = %d',1);
%
nip = fscanf(fid, '\nNUMBER_OF_INTEGRATION_POINTS = %d',1);
%
ngaus = fscanf(fid, '\nNUMBER_OF_GAUSS_POINTS = %d',1);
%
ntype = fscanf(fid, '\nANALYSIS_TYPE = %d', 1);
%
% ntype = 1  ->  Plane stress analysis
% ntype = 2  ->  Plane strain analysis
% ntype = 3  ->  Axisymmetric analysis
%
% Total number of elements in the mesh
nelem = fscanf(fid, '\nELEMENTS = %d', 1);
%
% Table of connectivities
%
lnods = fscanf(fid, '\n%d', [1+nnode,nelem]);
lnods = lnods';sortrows(lnods,1);
lnods = lnods(:,2:nnode+1);
%...create table of element global degrees of freedom
eldofX = lnods*2-1; eldofY = lnods*2;
if nnode == 3
    eldofs = [ eldofX(:,1) eldofY(:,1) eldofX(:,2) eldofY(:,2) eldofX(:,3) eldofY(:,3) ];
elseif nnode == 4 
    eldofs = [ eldofX(:,1) eldofY(:,1) eldofX(:,2) eldofY(:,2) ...
               eldofX(:,3) eldofY(:,3) eldofX(:,4) eldofY(:,4) ];
elseif nnode == 8
    eldofs = [ eldofX(:,1) eldofY(:,1) eldofX(:,2) eldofY(:,2) ...
               eldofX(:,3) eldofY(:,3) eldofX(:,4) eldofY(:,4) ...
               eldofX(:,5) eldofY(:,5) eldofX(:,6) eldofY(:,6) ...
               eldofX(:,7) eldofY(:,7) eldofX(:,8) eldofY(:,8) ];
end
%
% Nodal coordinates
%
npoin = fscanf(fid, '\nNODE_COORDINATES = %d', 1);
coord = fscanf(fid, '\n%d %f %f', [3, npoin]);
coord = coord';sortrows(coord,1);nodnumbers=coord(:,1);
coord = coord(:,2:3);
%...create table of nodal degrees of freedom
nodofs = [ nodnumbers*2-1 nodnumbers*2 ];
%
% Material thickness (plane stress only)
%
if  ntype == 1
    thick = fscanf(fid, '\nTHICKNESS = %f', 1);
else
    thick = 0;
end
%
% Prescribed displacements
%
nnodefix = fscanf(fid,'\nNODES_WITH_PRESCRIBED_DISPLACEMENTS = %d',1);
ipresc   = fscanf(fid, '\n%d %d %f %f', [4, nnodefix]);
ipresc   = ipresc';

%...create tables of fixed dofs and corresponding prescribed values
ifdoffix = zeros(npoin*ndofn,1); icount = 0;
for inodefix = 1:nnodefix
    ipoin=ipresc(inodefix,1);
    dofX=nodofs(ipoin,1);dofY=nodofs(ipoin,2);
    if  ipresc(inodefix,2)==11
        ifdoffix(dofX)=1;
        ifdoffix(dofY)=1; 
        icount = icount+1;
        valuedoffix(icount)   = ipresc(inodefix,3);
        fixeddoftable(icount) = dofX;
        icount = icount+1;
        valuedoffix(icount)   = ipresc(inodefix,4);
        fixeddoftable(icount) = dofY;
    elseif  ipresc(inodefix,2)==1
        ifdoffix(dofY)=1;
        icount = icount+1;
        valuedoffix(icount)   = ipresc(inodefix,4);
        fixeddoftable(icount) = dofY;
    elseif  ipresc(inodefix,2)==10
        ifdoffix(dofX)=1;
        icount = icount+1;
        valuedoffix(icount)   = ipresc(inodefix,3);
        fixeddoftable(icount) = dofX;
    elseif  ipresc(inodefix,2)==0
    else
        error('Wrong displacement prescription code in data file')
    end
end
%...create table of free dofs by subtracting the set of fixed dofs
%   from the set of all dofs
ngdof=npoin*ndofn;alldoftable=zeros(ngdof,1);
for igdof=1:ngdof
    alldoftable(igdof)=igdof;
end
freedoftable = setxor(alldoftable,fixeddoftable);
%
% Material properties
%
matprop.young = fscanf(fid, '\nYOUNG_MODULUS = %f', 1);
matprop.poiss = fscanf(fid, '\nPOISSON_RATIO = %f', 1);
%
% Compute elasticity D matrix
%
young = matprop.young; poiss = matprop.poiss;
if  ntype == 1                   % plane stress
    
    const = young/(1-poiss^2);
    dmatx = const*[  1    poiss       0      ;...
                   poiss    1         0      ;...
                     0      0    (1-poiss)/2 ];
             
elseif  ntype == 2 | ntype == 3 % Plane strain, Axisymmetric
    
    const = young/((1+poiss)*(1-2*poiss));
    
    dmatx = const*[1-poiss   poiss          0        ;...
                    poiss   1-poiss         0        ;...
                      0        0       (1-2*poiss)/2 ];

    if ntype == 3 % add extra row for axisymmetric case

       dmatx = [    dmatx    const*[poiss  poiss  0]' ;...
                   const*[poiss  poiss   0   1-poiss] ];
    end
   
end
%
% Load vector (point loads only)
%
npload = fscanf(fid,'\nPOINT_LOADS = %d',1);
pload  = fscanf(fid, '\n%d %f %f', [3, npload]);
pload  = pload';
%...add point loads to global load vector
F1=zeros(ngdof,1);loadednodes=pload(:,1);loadeddofs=nodofs(loadednodes,:);
F1(loadeddofs) = pload(:,2:3);
% elements in part 1;
p1 = [1  2  3  4 13 14 15 16 25 26 27 28 37 38 39 40];
p2 = [5  6  7  8 17 18 19 20 29 30 31 32 41 42 43 44];
p3 = [9 10 11 12 21 22 23 24 33 34 35 36 45 46 47 48];
l1 = unique(eldofs(p1,:)); K1 = zeros(length(l1),length(l1));
l2 = unique(eldofs(p2,:)); K2 = zeros(length(l2),length(l2));
l3 = unique(eldofs(p3,:)); K3 = zeros(length(l3),length(l3));
p4 = [1  2  3  4 13 14 15 16 25 26 27 28 37 38 39 40 5 6 17 18 29 30 41 42];
p5 = [9 10 11 12 21 22 23 24 33 34 35 36 45 46 47 48 7 8 19 20 31 32 43 44];
l4 = unique(eldofs(p4,:)); K4 = zeros(length(l4),length(l4));
l5 = unique(eldofs(p5,:)); K5 = zeros(length(l5),length(l5));
pa = [  1  2  3  4 13 14 15 16 ]; 
pb = [  5  6  7  8 17 18 19 20 ];
pc = [  9 10 11 12 21 22 23 24 ];
pd = [ 25 26 27 28 37 38 39 40 ];
pe = [ 29 30 31 32 41 42 43 44 ];
pf = [ 33 34 35 36 45 46 47 48 ];
la = unique(eldofs(pa,:)); Ka = zeros(length(la),length(la));
lb = unique(eldofs(pb,:)); Kb = zeros(length(lb),length(lb));
lc = unique(eldofs(pc,:)); Kc = zeros(length(lc),length(lc));
ld = unique(eldofs(pd,:)); Kd = zeros(length(ld),length(ld));
le = unique(eldofs(pe,:)); Ke = zeros(length(le),length(le));
lf = unique(eldofs(pf,:)); Kf = zeros(length(lf),length(lf));
%
% Solution phase
% ==============
%
% Compute global stiffness and assemble system of equations
%
K = zeros(npoin*ndofn,npoin*ndofn);
for ielem = 1:nelem
%...compute element stiffness
    if nnode == 3
        ke = stiffT3(coord,ielem,lnods,ntype,thick,dmatx);
    elseif nnode == 4 || nnode == 8
        ke = stiffQ(coord,ielem,lnods,ntype,thick,dmatx,nip,ngaus);
    end
%...add element contribution to global stiffness
    gpos = eldofs(ielem,:);
    K(gpos,gpos) = K(gpos,gpos) + ke;
    if ismember(ielem,p1)
        [~,target_dof] = ismember(eldofs(ielem,:),l1);
        K1(target_dof,target_dof) = K1(target_dof,target_dof) + ke;
    elseif ismember(ielem,p2)
        [~,target_dof] = ismember(eldofs(ielem,:),l2);
        K2(target_dof,target_dof) = K2(target_dof,target_dof) + ke;
    elseif ismember(ielem,p3)
        [~,target_dof] = ismember(eldofs(ielem,:),l3);
        K3(target_dof,target_dof) = K3(target_dof,target_dof) + ke;
    end
    if ismember(ielem,p4)
        [~,target_dof] = ismember(eldofs(ielem,:),l4);
        K4(target_dof,target_dof) = K4(target_dof,target_dof) + ke;
    elseif ismember(ielem,p5)
        [~,target_dof] = ismember(eldofs(ielem,:),l5);
        K5(target_dof,target_dof) = K5(target_dof,target_dof) + ke;
    end
    if ismember(ielem,pa)
        [~,target_dof] = ismember(eldofs(ielem,:),la);
        Ka(target_dof,target_dof) = Ka(target_dof,target_dof) + ke;
    elseif ismember(ielem,pb)
        [~,target_dof] = ismember(eldofs(ielem,:),lb);
        Kb(target_dof,target_dof) = Kb(target_dof,target_dof) + ke;
    elseif ismember(ielem,pc)
        [~,target_dof] = ismember(eldofs(ielem,:),lc);
        Kc(target_dof,target_dof) = Kc(target_dof,target_dof) + ke;
    elseif ismember(ielem,pd)
        [~,target_dof] = ismember(eldofs(ielem,:),ld);
        Kd(target_dof,target_dof) = Kd(target_dof,target_dof) + ke;
    elseif ismember(ielem,pe)
        [~,target_dof] = ismember(eldofs(ielem,:),le);
        Ke(target_dof,target_dof) = Ke(target_dof,target_dof) + ke;
    elseif ismember(ielem,pf)
        [~,target_dof] = ismember(eldofs(ielem,:),lf);
        Kf(target_dof,target_dof) = Kf(target_dof,target_dof) + ke;
    end
end
K;
%...assemble reduced stifness matrix and load vector
%   (remove prescibed d.o.f's)
Kstar = K(freedoftable,freedoftable);
Fstar = F1(freedoftable);
%...add contributions from prescribed displacements to
%   right hand side (Fstar)
Fstar=Fstar-K(freedoftable,fixeddoftable)*valuedoffix';
%
% Solve system of equations for unknown displacement
% vector Ustar
%
Ustar = Kstar \ Fstar;
%...assemble full global displacement vector
U=zeros(ngdof,1);U(freedoftable)=Ustar;
U(fixeddoftable)=valuedoffix;
% [U(1:2:length(U)),U(2:2:length(U))]
%
% Compute nodal reactions
%
R=zeros(ngdof,1);
R(fixeddoftable)=K(fixeddoftable,:)*U-F1(fixeddoftable);
% [R(1:2:length(R)),R(2:2:length(R))]


dummyU = U;
% common_ab = intersect(la,lb);
% [~,targeta] = ismember(common_ab,la);
% [~,targetb] = ismember(common_ab,lb);
% [Ka(targeta,:)*U(la) Kb(targetb,:)*U(lb)]
% dummyU(common_ab(1:2:length(common_ab))) = [ 0.18; 0.22; 0.18 ];
% [Ka(targeta,:)*U(la) Kb(targetb,:)*U(lb)]
% F_c_a = Ka(targeta,:)*dummyU(la);
% F_c_b = Ka(targetb,:)*dummyU(lb);
% [F_c_a  F_c_b]
% % for element a
% dof_a_i = [17;18];  dof_a_b = setdiff(la,dof_a_i); 
% [~,ta] = ismember(dof_a_b,la);
% Fa = Ka*dummyU(la); 
% Fa_b = Fa(ta);
% dof_b_i = [21;22];  dof_b_b = setdiff(lb,dof_b_i);
% [~,tb] = ismember(dof_b_b,lb);
% Fb = Kb*dummyU(lb); 
% Fb_b = Fb(tb);
% Fext = (Ka(targeta,:)*dummyU(la)+Kb(targetb,:)*dummyU(lb))/2;
% zerodof = find(abs(Fext) < 1e-10);
% nonzore = setdiff([1:length(Fext)],zerodof);
% ryvrsa = targeta(zerodof);
% ryvrsb = targetb(zerodof);
% F_c_a(zerodof) = zeros(length(ryvrsa),1);
% F_c_b(zerodof) = zeros(length(ryvrsb),1);
% resa = setdiff(targeta,ryvrsa);
% resb = setdiff(targetb,ryvrsb);
% F_c_a(nonzore) = F_c_a(nonzore) - Fext(nonzore);
% F_c_b(nonzore) = F_c_b(nonzore) - Fext(nonzore);
% [F_c_a F_c_b]
common_12 = intersect(l1,l2);
[~,target1] = ismember(common_12,l1);
[~,target2] = ismember(common_12,l2);
[K1(target1,:)*U(l1) K2(target2,:)*U(l2)]
dummyU(common_12(1:2:length(common_12))) = [ 0.28; 0.32; 0.28; 0.32; 0.28 ];
[K1(target1,:)*dummyU(l1) K2(target2,:)*dummyU(l2)]
F_c_1 = K1(target1,:)*dummyU(l1);
F_c_2 = K2(target2,:)*dummyU(l2);
[F_c_1  F_c_2]
% for element a
dof_1_i = [17;18];  dof_1_b = setdiff(l1,dof_1_i); 
[~,t1] = ismember(dof_1_b,l1);
F1 = K1*dummyU(l1); 
dof_2_i = [21;22];  dof_2_b = setdiff(l2,dof_2_i);
[~,tb] = ismember(dof_2_b,l2);
F2 = K2*dummyU(l2); 
Fext = (K1(target1,:)*dummyU(l1)+K2(target2,:)*dummyU(l2))/2;
zerodof = find(abs(Fext) < 1e-10);
nonzore = setdiff([1:length(Fext)],zerodof);
ryvrs1 = target1(zerodof);
ryvrs2 = target2(zerodof);
F_c_1(zerodof) = zeros(length(ryvrs1),1);
F_c_2(zerodof) = zeros(length(ryvrs2),1);
F_c_1(nonzore) = F_c_1(nonzore) - Fext(nonzore);
F_c_2(nonzore) = F_c_2(nonzore) - Fext(nonzore);
[F_c_1 F_c_2]

%
% Post-processing phase
% =====================
%
%...split U and R into x and y components for plotting
Ux = zeros(ngdof/2,1);
Uy = zeros(ngdof/2,1);
Ry = zeros(ngdof/2,1);
Rx = zeros(ngdof/2,1);
for igxdof = 1:ngdof/2
    Ux(igxdof) = U(igxdof*2-1);
    Uy(igxdof) = U(igxdof*2);
    Rx(igxdof) = R(igxdof*2-1);
    Ry(igxdof) = R(igxdof*2);
end
%
% Plot original and deformed mesh (with amplified displacements)
%
Adjacency = sparse(npoin,npoin);
if nnode == 3
    matxposn12=lnods(:,1)+npoin*(lnods(:,2)-1);
    matxposn23=lnods(:,2)+npoin*(lnods(:,3)-1);
    matxposn31=lnods(:,3)+npoin*(lnods(:,1)-1);
    matxpos=[matxposn12;matxposn23;matxposn31];
elseif nnode == 4
    matxposn12=lnods(:,1)+npoin*(lnods(:,2)-1);
    matxposn23=lnods(:,2)+npoin*(lnods(:,3)-1);
    matxposn34=lnods(:,3)+npoin*(lnods(:,4)-1);
    matxposn41=lnods(:,4)+npoin*(lnods(:,1)-1);
    matxpos=[matxposn12;matxposn23;matxposn34;matxposn41];
end
Adjacency(matxpos)=1;
%...reads deformation amplification factor for plotting
amplfact = fscanf(fid,'\nPLOTTING_AMPLIFICATION_FACTOR = %d',1);
%...compute deformed coordinates (with amplification)
defcoord(:,1) = coord(:,1)+amplfact*Ux;
defcoord(:,2) = coord(:,2)+amplfact*Uy;
%...and finally plots undeformed and deformed meshes
subplot(1,2,1);
gplot(Adjacency,coord,'-*b');hold on;
gplot(Adjacency,defcoord,'-r');hold off;
%
% Plot nodal reactions on top of undeformed mesh
%
subplot(1,2,2);
gplot(Adjacency,coord,'-b');hold on;
quiver(coord(:,1),coord(:,2),Rx,Ry,'r');hold off;
%
% Close file(s) and clear all variables before terminating program
% ----------------------------------------------------------------
%
status = fclose(fid);
clear;
