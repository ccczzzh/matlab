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
F=zeros(ngdof,1);loadednodes=pload(:,1);loadeddofs=nodofs(loadednodes,:);
F(loadeddofs) = pload(:,2:3);
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
end
K;
%...assemble reduced stifness matrix and load vector
%   (remove prescibed d.o.f's)
Kstar = K(freedoftable,freedoftable);
Fstar = F(freedoftable);
%...add contributions from prescribed displacements to
%   right hand side (Fstar)
Fstar=Fstar-K(freedoftable,fixeddoftable)*valuedoffix';
%
% Solve system of equations for unknown displacement
% vector Ustar
%
Ustar = Kstar \ Fstar;
%...assemble full global displacement vector
U=zeros(ngdof,1);U(freedoftable)=Ustar;U(fixeddoftable)=valuedoffix
%
% Compute nodal reactions
%
R=zeros(ngdof,1);
R(fixeddoftable)=K(fixeddoftable,:)*U-F(fixeddoftable);
R;
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
%clear;
