function [area,ke,intbmatx,lnods_trans,accbmatx] = stiffness(dmatx,domain,ngaus,coord,ielem,lnods,ntype,thick,nnode,ndim,ielgroup)
%************************************************************************
% Depending on the FE Method and Gauss Quadrature, computes the stiffness
% for all kinds of elements(T3, T6, Q4, Q8, Q9) in plane stress, plane
% strain and axisymmetric elasticity conditions. Also according to the
% multi-scale modeling, accumulation of the B matrix is done in this part.
%
% INPUT ARGUMENTS
% dmatx      - D matrix.
% domain     - Domain types:
%              'QUA' for quadrilateral domains.
%              'TRI' for triangular domains.
% ngaus      - Number of guass points in one elements.
% coord      - Coordinates of each node.
% ielem      - The current looping element (according to 'microcell')
% lnods      - Local nodes numbers that consists elements.
% ntype      - Analysis type
%                  ntype = 1  ->  Plane stress analysis
%                  ntype = 2  ->  Plane strain analysis.
% thick      - Thickness of the RVE.
% nnode      - Number of nodes concluded in one element.
% ndim       - Number of spatial dimensions.
%
% OUTPUT
% area       - Area of the current looping element.
% ke         - The stiffness matrix for the current looping element.
% intbmatx   - integral of B matrix in one element.
%
%***********************************************************************
ndofn =2 ; nevab = nnode*ndofn; area = 0; intbmatx=zeros(3,nevab);
lnods_trans = zeros(1,nnode); accbmatx = zeros(3,nevab);
% Loop the gauss points
ke = zeros(nevab,nevab);
if nnode == 3; % distinguish the element type through nnode
    x = zeros(nnode,1);
    y = zeros(nnode,1);
    for inode = 1:nnode
        x(inode)=coord(lnods(ielem,inode),1);
        y(inode)=coord(lnods(ielem,inode),2);
    end
    detj=(x(1)-x(3))*(y(2)-y(3))-(y(1)-y(3))*(x(2)-x(3));
    % Sometimes the elements could be built clockwise or anti-clockwise,
    % which could cause the detj < 0, while detj < 0, just reverse the
    % lnods.
    if detj < 0
        [lnods_trans(1,:)] = trans(lnods, nnode, ielem);
        lnods(ielem,:) = lnods_trans;
        for inode = 1:nnode
            x(inode)=coord(lnods_trans(1,inode),1);
            y(inode)=coord(lnods_trans(1,inode),2);
        end
    end
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
    if ntype == 1
        dvolu=thick*area;
    elseif ntype == 2
        dvolu=area;
    elseif ntype == 3
        twopi=8*atan(1);
        dvolu=twopi*rc*area;
    end
    ke = bmatx'*dmatx*bmatx*dvolu;
    intbmatx= intbmatx + bmatx*dvolu;
elseif nnode == 4 || nnode ==6 || nnode ==8 || nnode == 9;
    [posgp,weigp] = gaus2d(domain,ngaus);
    for igaus = 1:ngaus;
        s=posgp(1,igaus);
        t=posgp(2,igaus);
        [deriv] = SFandD(s,t,nnode);
        %
        % assemble the local node number and their corresponding coordinates
        %
        eloc = zeros(nnode,2);
        for in = 1:nnode;
            eloc(in,:)= coord(lnods(ielem,in),:);
        end
        % calculate the jacobian matrix
        jac=deriv*eloc;
        detJ=det(jac);
        if detJ<=0
            % The reason that the determinent is negative is the lnods is not
            % coupled with the derivative, in this code, only one other option
            % is tried
            [lnods_trans(1,:)] = trans(lnods, nnode, ielem);
            lnods(ielem,:) = lnods_trans;
            for in = 1:nnode;
                eloc(in,:)= coord(lnods_trans(1,in),:);
            end
            jac=deriv*eloc;
            detJ=det(jac);
            if detJ<0;
                error('Error: Negative Jacobian, the Mesh is not well organised.');
            end
        end
        jaci=inv(jac);
        % calculate the cartisian product
        cartd=jaci*deriv;
        % build up the b matrix for the looping element
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
        % Calculate area
        iarea= detJ*weigp(igaus);%
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
        ke = ke+bmatx'*dmatx(:,:,ielgroup)*bmatx*dvolu;
        area = area + iarea;
        %accumulate the bmaxt
        intbmatx = intbmatx + bmatx*dvolu;
        accbmatx = accbmatx + bmatx;
    end
end
lnods_trans = lnods(ielem,:);
end

