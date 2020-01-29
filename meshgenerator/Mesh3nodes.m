function [coords nT nNodes ]=Mesh3nodes(Lx,Ly,Nx,Ny)
%   This function generates triangular mesh for a rectangular
%   shape structure for finite element analysis
%   [coords nT nNodes ]=femTriangularMeshGenerator(Lx,Ly,Nx,NE)
%   coords  =   x and y coordinates of each element nodes
%   nT      =   nodal connectivity
%   nNodes  =   Number of nodes
%   Lx      =   width of the rectangular structure
%   Ly      =   Height of the rectangular structure
%   Nx      =   Number of divisions on x- axis
%   NE      =   Number of elements
%



NE = 2*Ny*Nx;
nNodes = (Nx+1)*(Ny+1);

m=1;
j=1:Nx;
k=linspace(Nx*2,NE,Ny);

s = Ny;

for i=1:s
    nT(m:2:k(i),1)= j;  %node 1 of 1st element
    nT(m+1:2:k(i),1)= j;  %node 1 of 2nd element
    nT(m:2:k(i),2)=j+1; %%node 2 of 1st element
    nT(m+1:2:k(i),2)=j+Nx+2;%node 2 of 2nd element
    nT(m:2:k(i),3)= j+Nx+2;  %node 3 of 1st element
    nT(m+1:2:k(i),3)=j+1+Nx;    %%node 3 of 2nd element
    
    m=k(i)+1;
    j=j+Nx+1;
end
%%%%%%%%%%%%%%%COORDINATES GENERATION%%%%%%%%%%%%%%%%%%%

ax=linspace(0,Lx,Nx+1); %%%x coordinates
by=linspace(0,Ly,Ny+1); %%%y coordinates
X1=[];
Y1=[];
for i1=1:Ny+1
    %    General Nodal Coordinates layer by layer
    by1(1:Nx+1)=by(i1);
    X1=[X1 ax];
    Y1=[Y1 by1];
end
j=1:3;

%each element coordinates
for n=1:NE
    X(j,1) = X1(nT(n,:));
    Y(j,1)=Y1(nT(n,:));
    j=j+3;
end
coords=[X Y];   %x and y coordinates