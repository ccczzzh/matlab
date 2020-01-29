%**********************************************************************
% Program for Dynamic Finite Element Analysis of Linear Elastic Solids
% in Plane Stress, Plane Strain and Axisymmetric Conditions with
% 3-noded Triangular Elements.
% Uses the Newmark Method for time integration
%
% HISTORY
% EAdSN  March 2015:  Initial coding, starting from "elast2D_v2"
%***********************************************************************
clear;cla;hold off;
% Set some basic variables
% ========================
% Number of degrees of freedom per node
ndofn=2;
% Number of nodes of the element
nnode=3;
%
% Data input phase. Read information from data file
% =================================================
%
filename = input('Input data file name -------> ','s')
fid = fopen(filename, 'r');

title = fscanf(fid, 'TITLE = %s',1)
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
lnods = lnods';
lnods = lnods(:,2:nnode+1);
%...create table of element global degrees of freedom
eldofX = lnods*2-1; eldofY = lnods*2;
eldofs = [ eldofX(:,1) eldofY(:,1) eldofX(:,2) eldofY(:,2) eldofX(:,3) eldofY(:,3) ];
%
% Nodal coordinates and velocity initial conditions
%
npoin = fscanf(fid, '\nNODE_COORDINATES_AND_INITIAL_VELOCITIES = %d', 1);
coord = fscanf(fid, '\n%d %f %f %f %f', [5, npoin]);
coord = coord'; 
nodnumbers=coord(:,1);
vinit = coord(:,4:5); coord = coord(:,2:3);
%...create table of nodal degrees of freedom
nodofs = [ nodnumbers*2-1 nodnumbers*2 ];
%...arrange vector of initial nodal velocities in usual manner
Vinit(nodofs(:,1))=vinit(:,1);Vinit(nodofs(:,2))=vinit(:,2);Vinit=Vinit';
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
%
%...create tables of fixed dofs and corresponding prescribed values
%   Note: actual prescriptions will be the values prescribed here times
%         the corresponding load factor (a prescribed function of time)
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
ngdof=npoin*ndofn;alldoftable=1:ngdof;
freedoftable = setxor(alldoftable,fixeddoftable);
%
% Function defining the load factor for prescribed displacements as a
% function of time (piecw-wice linear function)
ndlfpoints = fscanf(fid,'\nLOAD_FACTOR_FUNCTION_FOR_PRESCRIBED_DISPLACEMENTS = %d',1);
dlffunc    = fscanf(fid, '\n%f %f', [2, ndlfpoints]);
%
% Material properties
%
matprop.density = fscanf(fid, '\nMASS_DENSITY = %f', 1);
matprop.young = fscanf(fid, '\nYOUNG_MODULUS = %f', 1);
matprop.poiss = fscanf(fid, '\nPOISSON_RATIO = %f', 1);
%
% Load vector (point loads only)
% Note: actual applied load will be the load prescribed here times the
%       relevant load factor (a prescribed function of time)
%
npload = fscanf(fid,'\nPOINT_LOADS = %d',1);
pload  = fscanf(fid, '\n%d %f %f', [3, npload]);
pload  = pload';
%...add point loads to global load vector
F=zeros(ngdof,1);loadednodes=pload(:,1);loadeddofs=nodofs(loadednodes,:);
F(loadeddofs) = pload(:,2:3);
%
% Function defining the load factor for loads as a function of time
% (piece-wise linear function)
nllfpoints = fscanf(fid,'\nLOAD_FACTOR_FUNCTION_FOR_PRESCRIBED_LOADS = %d',1);
llffunc    = fscanf(fid, '\n%f %f', [2, nllfpoints]);
%
%
% Analysis time interval
%
tfinal = fscanf(fid,'\nFINAL_TIME = %f',1);
% time step
dtime  = fscanf(fid,'\nTIME_STEP = %f',1);
%
% Newmark Method parameters for time integration:
%       0 <= gamma <= 1
%       0 <= beta  <= 1/2
%
%       Particular cases
%       - Constant average acceleration: gamma=1/2; beta=1/4
%       - Linear acceleration:  gamma=1/2; beta=1/6 (conditionally stable)
%
gamma = fscanf(fid,'\nNEWMARK_GAMMA_PARAMETER = %f',1);
beta  = fscanf(fid,'\nNEWMARK_BETA_PARAMETER = %f',1);
%
%
amplfac  = fscanf(fid,'\nPLOTTING_AMPLIFICATION_FACTOR = %f',1);
movdtime  = fscanf(fid,'\nTIME_INTERVAL_BETWEEN_MOVIE_FRAMES = %f',1);
%
% ==============
% Solution phase
% ==============
%
% Assemble global stiffness and mass matrices
%
K = zeros(npoin*ndofn,npoin*ndofn);M = zeros(npoin*ndofn,npoin*ndofn);
for ielem = 1:nelem
%...compute element stiffness and mass matrices
    [ke,me] = T3stiffmass(coord,ielem,lnods,matprop,ntype,thick);  
%...add element stiffness contribution to global stiffness
    gpos = eldofs(ielem,:);
    K(gpos,gpos) = K(gpos,gpos) + ke;
%...add element mass contribution to global (lumped) mass matrix
    M(gpos,gpos) = M(gpos,gpos) + me;
end
K;
M;
Cm = 76.6972;
Ck = 8.0835e-7;
C  = Cm*M + Ck*K;
%
% Time steping using Newmark time integration algorithm
% -----------------------------------------------------
%
% initialise time, nodal displacements, velocities and accelerations
time = 0;
V=Vinit;
dlfactor=plfunc(time,ndlfpoints,dlffunc);
U(fixeddoftable)=dlfactor*valuedoffix;
U(freedoftable)=0; % initial displacement of free dofs assumed zero
U=U';
A(fixeddoftable)=0;A=A'; % TEMPORARY CODE!!! needs fixing - assumes zero initial acceleration of prescribed dofs
%
% Solve for possible non-zero initial acceleration at free dofs
% (i.e. ensures the body is in dynamic equililibrium at t0)
%
llfactor=plfunc(time,nllfpoints,llffunc);
ForcingTerm = llfactor*F-K*U;
ForcingTerm = ForcingTerm(freedoftable);
ForcingTerm = ForcingTerm - M(freedoftable,fixeddoftable)*A(fixeddoftable);
%...solve for initial acceleration of free dofs
Mstar=M(freedoftable,freedoftable);
Astar = Mstar \ ForcingTerm;
A(freedoftable)=Astar;A(fixeddoftable)=0; % TEMPORARY CODE!!! needs fixing
%
% Proceed to subsequent time steps
%
ii=0;pltframe=0;accmdtime=0;
while time<=tfinal
    % update time, load factors and load vector
    time=time+dtime;accmdtime=accmdtime+dtime;
    llfactor=plfunc(time,nllfpoints,llffunc);
    Fnew=llfactor*F;
    dlfactor=plfunc(time,ndlfpoints,dlffunc);
    % re-set displacement, velocity and acceleration vectors
    Aold=A;Vold=V;Uold=U;
    % Solve dynamic equilibrium FE equations for current acceleration
    %...assemble system matrix and right-hand-side
    SysMtx=M+C*dtime*gamma+dtime^2*beta*K;
    SysRhs=Fnew-K*(Uold+dtime*Vold+dtime^2*(1/2-beta)*Aold)-C*(Vold+dtime*(1-gamma)*Aold);
    %...add contribution of current time acceleration of prescribed dofs to
    %   right-hand-side
    SysRhs=SysRhs(freedoftable); % remove prescribed dofs
    valueaccpresc=zeros(length(fixeddoftable),1); % TEMPORARY CODE!!! needs fixing
    SysRhs=SysRhs-SysMtx(freedoftable,fixeddoftable)*valueaccpresc;
    SysMtx=SysMtx(freedoftable,freedoftable); % remove prescribed dofs
    %...obtain current time acceleration of free dofs
    Astar = SysMtx \ SysRhs;
    % Update free velocities
    Vstar = Vold(freedoftable)+dtime*((1-gamma)*Aold(freedoftable)+gamma*Astar);
    % Update free displacements
    Ustar = Uold(freedoftable)+dtime*Vold(freedoftable)+dtime^2/2*((1-2*beta)*Aold(freedoftable)+beta*Astar);
    %...add prescribed dofs (assemble full vectors)
    A(freedoftable)=Astar;
    A(fixeddoftable)=0; % TEMPORARY CODE!!! needs fixing
    V(freedoftable)=Vstar;ddlfactor=dplfunc(time,ndlfpoints,dlffunc);
    V(fixeddoftable)=ddlfactor*valuedoffix;
    U(freedoftable)=Ustar;U(fixeddoftable)=dlfactor*valuedoffix;
    ii=ii+1;UU(ii)=U(5);tt(ii)=time;
    Ux=U(1:2:npoin*ndofn);Uy=U(2:2:npoin*ndofn);
    % plot frame at or near the pre-specified movie frame time station
    if accmdtime >= movdtime
        subplot(1,2,1);hold off;
        patch('Faces',lnods,'Vertices',coord+amplfac*[Ux Uy],'FaceColor','W','EdgeColor','k');
        axis equal;
        MM(ii)=getframe;
        accmdtime=0; % re-set accumulated time interval
    end
end
%movie(MM,1);
%
% =====================
% Post-processing phase
% =====================
%
subplot(1,2,2);plot(tt,UU);
%
%
hold off;
%
% Close file(s) and clear all variables before terminating program
% ----------------------------------------------------------------
%
status = fclose(fid);
clear;
