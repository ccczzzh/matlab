%
% Computes the stiffness for 3-noded triangular elements
% in plane stress, plane strain and axisymmetric elasticity
%
function [ke] = stiffT3(coord,ielem,lnods,ntype,thick,dmatx)
%
%
nnode = 3; ndofn =2 ; nevab = nnode*ndofn;
% Retrieve element nodal coordinates
x = zeros(nnode,1);
y = zeros(nnode,1);
for inode = 1:nnode
    x(inode)=coord(lnods(ielem,inode),1);
    y(inode)=coord(lnods(ielem,inode),2);
end
% Compute element area
detj=(x(1)-x(3))*(y(2)-y(3))-(y(1)-y(3))*(x(2)-x(3));
area=detj/2;
%
% Assemble strain-displacement B matrix
%
%... part of B common to plane stress,
%    plane strain and axisymmetric
bmatx=[y(2)-y(3)      0      y(3)-y(1)      0      y(1)-y(2)      0    ;...
           0      x(3)-x(2)      0      x(1)-x(3)     0       x(2)-x(1);...
       x(3)-x(2)  y(2)-y(3)  x(1)-x(3)  y(3)-y(1)  x(2)-x(1)  y(1)-y(2)]...
      /detj;
%...add extra row for axisymmetric case
if ntype == 3  
    rc = (x(1)+x(2)+x(3))/3; % centroid radial coordinate
    a  = 1/(3*rc);
    bmatx = [          bmatx        ;...
              a   0   a   0   a   0 ];
end
%
% Compute element stiffness
%
%...differential volume
if ntype == 1
   dvolu=thick*area;
elseif ntype == 2
   dvolu=area;
elseif ntype == 3
   twopi=8*atan(1);
   dvolu=twopi*rc*area;
end
%...element stiffness
ke = bmatx'*dmatx*bmatx*dvolu;
