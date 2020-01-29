%
% Computes the stiffness for 4-noded quadratic elements
% in plane stress, plane strain and axisymmetric elasticity
%
function [keg] = stiffQG(coord,ielem,lnods,ntype,thick,dmatx,nip,ngaus)
%
%
nnode = 4; ndofn =2; nevab = nnode*ndofn; area = 0; domain = 'QUA';
ndim = 2; 
keg = zeros(nevab,nevab);
% for gauss integration
[posgp,weigp] = gaus2d(domain,ngaus);
for igaus = 1:ngaus
    s=posgp(1,igaus);
    t=posgp(2,igaus);
    [deriv] = SFandD(s,t,nnode);
    %
    % assemble the local node number and their corresponding coordinates
    %
    eloc = zeros(nnode,2);
    for in = 1:nnode
        eloc(in,:)= coord(lnods(ielem,in),:);
    end
    % calculate the jacobian matrix
    jac=deriv*eloc;
    detJ=det(jac);
    jaci=inv(jac);
    % calculate the cartisian product
    cartd=jaci*deriv;
    bmatx=zeros(3,nnode*ndofn);
    for im=1:nnode
        i=(im-1)*ndim+1 ; j=i+1  ;
        bmatx(1,i)=cartd(1,im);
        bmatx(1,j)=0.0;
        bmatx(2,i)=0.0;
        bmatx(2,j)=cartd(2,im);
        bmatx(3,i)=cartd(2,im);
        bmatx(3,j)=cartd(1,im);
    end
%     % Calculate area
%     iarea= detJ*weigp(igaus);
    % Calculate area
    iarea= detJ*weigp(igaus);
    % Compute element stiffness
    %
    %...differential volume
    if ntype == 1
        dvolu=thick*iarea;
    elseif ntype == 2
        dvolu=iarea;
    elseif ntype == 3
        twopi=8*atan(1);
        dvolu=twopi*rc*iarea;
    end
    %...element stiffness, use gauss quadrature integration, so its looping
    % aboud the gauss points
    keg = keg+bmatx'*dmatx*bmatx*dvolu;
    area = area + iarea;
    %accumulate the bmaxt
end

