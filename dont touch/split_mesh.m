function [bound,cell_volume,dofs,nodes] = split_mesh(coord,lnods_final,nelem,ndofn,npoin,option,nnode)
%***********************************************************************
% function split_mesh
%
% DESCRIPTION
% Splits the nodes and degrees of freedom of a given rectangular micro-
% cell mesh into sets of nodes and degrees of freedom associated with
% interior, boundary and different portions of the boundary.
% It also checks for one-to-one correspondance between nodes on opposing
% sides of the micro-cell rectangle, which is required in the present
% computational treatment of models under the assumption of periodic
% boundary displacement fluctuations with anti-periodic tractions.
% Also creates the structure "bound" containing information
% needed for the uniform boundary traction case.
%
% Still, the RVE is concidered as rectangulars, which has straight
% boundaries and boundaries are either horizontal or vertical.
%
% INPUT ARGUMENTS
% coord                - coordinates of nodes
% lnods                - local node numbers that consists elements
% nelem                - number of elements
% ndofn                - number of degree of freedom
% nopin                - total number of nodes
% option               - RVE kinematics assumption option
% nnode                - number of nodes contained one element
%
% OUTPUT
% bound.R              - the calculation of -inv(Cd)*Cf;
%
%***********************************************************************
% Find maximum and minimum coordinates (coordinates of sides of
% rectangular micro-cell), dimensions and volume of micro-cell
maxcoord    = max(coord);
mincoord    = min(coord);
cell_length = maxcoord - mincoord;
cell_volume = prod(cell_length);
%
% Find nodes on micro-cell boundary and create arrays (sets) of boundary
% and interior node numbers
% ----------------------------------------------------------------------
%
toler = 10^(-8)*cell_length;
maxmt = maxcoord-toler;
maxpt = maxcoord+toler;
minmt = mincoord-toler;
minpt = mincoord+toler;
node_corner = zeros(4,1);
node_interior = 0; node_edge_left   = 0; node_edge_right = 0;
node_edge_top = 0; node_edge_bottom = 0;
for ipoin = 1:npoin
    bmin = coord(ipoin,:)>minmt & coord(ipoin,:)<minpt;
    bmax = coord(ipoin,:)>maxmt & coord(ipoin,:)<maxpt;
    if bmin == [0 0] & bmax == [0 0]
        node_interior = [node_interior; ipoin];
    elseif bmin == [1 1] & bmax == [0 0]
        node_corner(1) = ipoin;
    elseif bmin == [0 1] & bmax == [1 0]
        node_corner(2) = ipoin;
    elseif bmin == [0 0] & bmax == [1 1]
        node_corner(3) = ipoin;
    elseif bmin == [1 0] & bmax == [0 1]
        node_corner(4) = ipoin;
    elseif bmin == [1 0] & bmax == [0 0]
        node_edge_left   = [node_edge_left;   ipoin];
    elseif bmin == [0 0] & bmax == [1 0]
        node_edge_right  = [node_edge_right;  ipoin];
    elseif bmin == [0 1] & bmax == [0 0]
        node_edge_bottom = [node_edge_bottom; ipoin];
    elseif bmin == [0 0] & bmax == [0 1]
        node_edge_top    = [node_edge_top;    ipoin];
    else
        error('Something is wrong with boundary nodes of a given micro-cell mesh')
    end
end
ni                = length(node_interior)-1;
node_interior     = node_interior(2:ni+1);
nl                = length(node_edge_left)-1;
node_edge_left    = node_edge_left(2:nl+1);
nr                = length(node_edge_right)-1;
node_edge_right   = node_edge_right(2:nr+1);
nb                = length(node_edge_bottom)-1;
node_edge_bottom  = node_edge_bottom(2:nb+1);
nt                = length(node_edge_top)-1;
node_edge_top     = node_edge_top(2:nt+1);
node_corner_bottom_left  = node_corner(1);
node_corner_bottom_right = node_corner(2);
node_corner_top_right    = node_corner(3);
node_corner_top_left     = node_corner(4);
%
% Perform various checks to ensure validity of micro-cell mesh and split
% nodes into interior, "plus", "minus" and corner sets and also create
% boundary node set
%
%... check that the given rectangular cell mesh has 4 corner nodes
corner_check = node_corner == zeros(4,1);
if any(corner_check)
    error('Wrong number of corner nodes detected in micro-cell')
end
if strcmp(option,'periodic')
    % Check one-to-one correspondance between boundary nodes when
    % periodic micro-cell boundary displacement fluctuation option is
    % selected. In this case, the sets of "plus" and "minus" boundary
    % nodes must contain matching pairs of nodes on opposing faces of the
    % micro-cell
    if nr ~= nl
        error('Number of nodes of right and left edges of micro-cell do not coincide')
    elseif nt ~= nb
        error('Number of nodes of top and bottom edges of micro-cell do not coincide')
    end
    % initialise array of right-left pairs
    node_right_left = zeros(nr,2);
    % check one-to-one correspondance and set right-left matching pairs
    % array
    tolery = toler(2);
    for i = 1:nr
        coorpt = coord(node_edge_right(i),2) + tolery;
        coormt = coord(node_edge_right(i),2) - tolery;
        count = 0;
        coordnj = zeros(1,2);
        for j = 1:nl
            coordnj = coord(node_edge_left(j),2);
            if coordnj < coorpt & coordnj > coormt
                count = count + 1;
                node_right_left(i,:) = [node_edge_right(i) node_edge_left(j)];
            end
        end
        if count ~= 1
            error('No one-to-one correspondance between micro-cell right and left boundary nodes')
        end
    end
    % check one-to-one correspondance and set top-bottom matching pairs
    % array
    tolerx = toler(1);
    for i = 1:nt
        coorpt  = coord(node_edge_top(i),1) + tolerx;
        coormt  = coord(node_edge_top(i),1) - tolerx;
        count   = 0;
        coordnj = zeros(1,2);
        for j = 1:nb
            coordnj = coord(node_edge_bottom(j),1);
            if coordnj < coorpt & coordnj > coormt
                count = count + 1;
                node_top_bottom(i,:) = [node_edge_top(i) node_edge_bottom(j)];
            end
        end
        if count ~= 1
            error('No one-to-one correspondance between micro-cell top and bottom boundary nodes')
        end
    end
    [lrl,crl] = size(node_right_left);
    [ltb,ctb] = size(node_top_bottom);
    nodes.p(1:lrl)          = node_right_left(:,1);
    nodes.p(lrl+1:lrl+ltb)  = node_top_bottom(:,1);
    nodes.m(1:lrl)         = node_right_left(:,2);
    nodes.m(lrl+1:lrl+ltb) = node_top_bottom(:,2);
else
    % For localisation options other than periodic boundary displacement
    % fluctuation, no one-to-one correspondance between nodes on opposing
    % faces is required. In this case, the sets of "plus" and "minus"
    % boundary nodes are created without any further checks and do not
    % in general comprise matching pairs of nodes
    nodes.p  = [ node_edge_right'  node_edge_top'    ];
    nodes.m = [ node_edge_left'   node_edge_bottom' ];
end
nodes.i   = node_interior';
nodes.c   = node_corner';
nodes.b   = [ nodes.p  nodes.m  nodes.c ];
nodes.all = [ nodes.i  nodes.b ];
%
% Create arrays (sets) of interior, boundary ("plus" and "minus") and
% corner nodes degrees of freedom
% -------------------------------------------------------------------
%
npoini = length(nodes.i);
for i = 1:npoini
    j = i*ndofn-1;
    dofs.i(j)   = nodes.i(i)*ndofn-1;
    dofs.i(j+1) = nodes.i(i)*ndofn;
end
npoinp = length(nodes.p);
for i = 1:npoinp
    j = i*ndofn-1;
    dofs.p(j)   = nodes.p(i)*ndofn-1;
    dofs.p(j+1) = nodes.p(i)*ndofn;
end
npoinm = length(nodes.m);
for i = 1:npoinm
    j = i*ndofn-1;
    dofs.m(j)   = nodes.m(i)*ndofn-1;
    dofs.m(j+1) = nodes.m(i)*ndofn;
end
npoinc = length(nodes.c);
for i = 1:npoinc
    j = i*ndofn-1;
    dofs.c(j)   = nodes.c(i)*ndofn-1;
    dofs.c(j+1) = nodes.c(i)*ndofn;
end
%
dofs.b   = [ dofs.p dofs.m dofs.c ];
dofs.all = [ dofs.i dofs.b ];
bound = [];
if strcmp(option,'traction')
    %
    % WARNING: This procedure is NOT FOOL PROOF !!
    % It assumes that element faces are linear (each boundary
    % face has 2 nodes only) and also that the mesh is
    % well defined (whatever that means ! ).
    %
    % Set dependent boundary dofs
    % (chosen as x and y-dof's of bottom right corner and
    %  x-dof of top right corner)
    dofs.d = dofs.c(3:5);
    % Set independent boundary dofs (assuming that top
    % left corner is fully constrained and top right
    % is constrained in y to remove rigid body mode)
    dofs.f = [ dofs.p dofs.m dofs.c(1:2) ];
    % Set prescribed dofs
    dofs.z = dofs.c(6:8);
    %
    % Create structure "bound" containing all relevant boundary
    % information. This information is used in the main program
    % to impose the (non-conventional) integral boundary displacement
    % constraint.
    %
    % ...detect elements with at least 2 nodes on the cell boundary
    %    (potentially, these have an edge on the cell boundary
    
    [m,n]       = size(lnods_final);
    bound.lnods = zeros(1,n); bound.nodpos  = zeros(1,n);
    bound.elno = 0          ; bound.nnod = 0;
    for ielem = 1:nelem
        flag      = 0;
        aux_pos   = zeros(1,n);
        aux_bound = zeros(1,n);
        for elnod = 1:n
            for ipoin = nodes.b
                if lnods_final(ielem,elnod) == ipoin
                    flag = flag + 1;
                    aux_bound(flag) = ipoin;
                    % identify node position
                    if ipoin == node_corner_bottom_left
                        aux_pos(flag) = 1;
                    elseif ipoin == node_corner_bottom_right
                        aux_pos(flag) = 3;
                    elseif ipoin == node_corner_top_right
                        aux_pos(flag) = 5;
                    elseif ipoin == node_corner_top_left
                        aux_pos(flag) = 7;
                    else
                        for i = node_edge_bottom'
                            if ipoin == i
                                aux_pos(flag) = 2;
                            end
                        end
                        for i = node_edge_right'
                            if ipoin == i
                                aux_pos(flag) = 4;
                            end
                        end
                        for i = node_edge_top'
                            if ipoin == i
                                aux_pos(flag) = 6;
                            end
                        end
                        for i = node_edge_left'
                            if ipoin == i
                                aux_pos(flag) = 8;
                            end
                        end
                    end
                end
            end
        end
        if flag >= 2
            bound.nodpos = [bound.nodpos; aux_pos];
            bound.lnods  = [bound.lnods; aux_bound];
            bound.elno   = [bound.elno ielem];
            bound.nnod   = [bound.nnod flag];
        end
    end
    nelb         = length(bound.elno);
    bound.nodpos = bound.nodpos(2:nelb,:);
    bound.lnods  = bound.lnods(2:nelb,:);
    bound.elno   = bound.elno(2:nelb);
    bound.nnod   = bound.nnod(2:nelb);
    
    % Below are the changed parts. From this point, the main logic is
    % similar compare to the original one, which is give a value to each 
    % boundary points, then boundary elements will assign a number to it.   
    % According to compare the assigned number with all the possible check
    % points, once coincide, we could know which boundary the element 
    % belongs, then carry on the following analsis; just more checking 
    % points and varying factors considered compare to the original
    % logic. 
    % The assumed nodal_value are shown below:
    %      15000          10         20000
    %          _ _ _ _ _ _ _ _ _ _ _
    %         |          |          |
    %         |          |          |
    %       1 |          |          | 100
    %         |_ _ _ _ _ |0 _ _ _ _ |
    %         |          |          |
    %         |          |          |
    %         |          |          |
    %         |_ _ _ _ _ |_ _ _ _ _ |
    %     10000        1000      15000
    % The reason that give such values is to make sure for each check point
    % there will be an unique number correspond to it.
    % To avoid round off error occurs due to using decimals, such 'large
    % scale' numbers are used
    % Different element types' situations are considered,however mixed
    % element type meshes are still not valid.
    
    % Below are the new check functions
    % The check points below are for nodes in an elements happens on only
    % one boundary.
    if nnode == 3;
        check.bottom = [11000 11100 31000 2000 2100 11010 22000 52000 2001 ...
            16000 16010 66000 16001 2010 25000 25010];
        check.right  = [15100 15110 65100 15101 200 210 50200 201 10200  ...
            1200 20100 20101 30100 21100 35000 35001];
        check.top    = [20010 20011 30010 21010 20 21 10020 1020 15020 120 ...
            50010 51010 65010 50110 70000 71000];
        check.left   = [50001 51001 65001 50110 2 1002 15002 10200 20002 ...
            12 10001 10101 30001 10011 60000 60100 ];
        check.all = [check.bottom check.right check.top check.left];
        
    elseif nnode == 6 ;
        check.bottom = [12000 12100 32000 12010 3000 3100 23000 3010 53000 ...
            3001 17000 17010 67000 17001 26000 26010];
        check.right  = [15200 15210 65200 15201 300 310 50300 301 10300 ...
            1300 20200 20201 30200 21200 35100 35101];
        check.top    = [20020 20021 3002 21020 30 31 10030 1030 15030 130 ...
            50020 51020 65020 50120 70010 71010];
        check.left   = [50002 51002 65002 50102 3 1003 25003 ...
            103 20003 13 10002 10102 30002 10012 60001 60101];
        check.all = [check.bottom check.right check.top check.left];
        
    elseif nnode == 4 ;
        check.bottom = [11000 2000 16000 25000];
        check.right  = [15100 200 20100 35000];
        check.top    = [20010 20 50010 70000];
        check.left   = [50001 2 10001 60000];
        check.all = [check.bottom check.right check.top check.left];
        
    elseif nnode == 8 || nnode == 9;
        check.bottom = [12000 3000 26000 17000];
        check.right  = [15200 300 20200 35100];
        check.top    = [20020 30 50020 70010];
        check.left   = [50002 3 10002 60001];
        check.all = [check.bottom check.right check.top check.left];
    end
    
    [p,q]   = size(bound.nodpos);
    aux1 = zeros(1,2); aux3 = 0;
    % unit normals
    nright = [1 ; 0]; nleft = -nright;
    ntop   = [0 ; 1]; nbottom = -ntop;
    bound.normals = [0 0];
    % Accumulate the point value each element assgned, then compare the
    % number with the check-points, we would know where the element is.
    for ib = 1:p;
        number = sum(bound.nodpos(ib,:)==1)*10000+sum(bound.nodpos(ib, :)==2)*1000+sum(bound.nodpos(ib, :)==3)*15000 ...
            +sum(bound.nodpos(ib, :)==4)*100+sum(bound.nodpos(ib, :)==5)*20000+sum(bound.nodpos(ib, :)==6)*10 ...
            +sum(bound.nodpos(ib, :)==7)*50000+sum(bound.nodpos(ib, :)==8)*1;
        
        zerovec = zeros(1,length(check.all));
        if (~isequal((number == check.all),zerovec))
            if nnode == 3 || nnode ==4;
                aux1 = [aux1 ; bound.lnods(ib,1:2)];
                aux3 = [aux3 ; bound.elno(ib)];
                if (~isequal((number == check.bottom),zeros(1,length(check.bottom))))
                    bound.normals = [bound.normals ; nbottom'];
                elseif (~isequal((number == check.right),zeros(1,length(check.right))))
                    bound.normals = [bound.normals ; nright'];
                elseif (~isequal((number == check.top),zeros(1,length(check.top))))
                    bound.normals = [bound.normals ; ntop'];
                elseif (~isequal((number == check.left),zeros(1,length(check.left))))
                    bound.normals = [bound.normals ; nleft'];
                end
            elseif nnode ==6 || nnode == 8 || nnode ==9;
                if (~isequal((number == check.bottom),zeros(1,length(check.bottom))))
                    if coord(bound.lnods(ib,1),1)-coord(bound.lnods(ib,2),1)>0;
                        % Due to an element's lnode might start from any
                        % corner in anti-clockwise sequence, if one element's
                        % lnods start from bottom-right corner, bound.lnods
                        % will get result of 1 2 3 will actually be B-R 
                        % corner then two nodes on the bottom boundary,
                        % therefore the correct aux will be 2-3, 3-1,
                        % otherwise the result will always be 1-2, 2-3.
                        % This is why need to do such judgement.
                        aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                        aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,1)];
                    elseif coord(bound.lnods(ib,1),1)-coord(bound.lnods(ib,2),1)<0;
                        aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                        aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    end
                    bound.normals = [bound.normals ; nbottom'];
                    bound.normals = [bound.normals ; nbottom'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif (~isequal((number == check.right),zeros(1,length(check.right))))
                    if coord(bound.lnods(ib,1),2)-coord(bound.lnods(ib,2),2)>0;
                        aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                        aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,1)];
                    elseif coord(bound.lnods(ib,1),2)-coord(bound.lnods(ib,2),2)<0;
                        aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                        aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    end
                    bound.normals = [bound.normals ; nright'];
                    bound.normals = [bound.normals ; nright'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif (~isequal((number == check.top),zeros(1,length(check.top))))
                    if coord(bound.lnods(ib,1),1)-coord(bound.lnods(ib,2),1)<0;
                        aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                        aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,1)];
                    elseif coord(bound.lnods(ib,1),1)-coord(bound.lnods(ib,2),1)>0;
                        aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                        aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    end
                    bound.normals = [bound.normals ; ntop'];
                    bound.normals = [bound.normals ; ntop'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif (~isequal((number == check.left),zeros(1,length(check.left))))
                    if coord(bound.lnods(ib,1),2)-coord(bound.lnods(ib,2),2)<0;
                        aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                        aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,1)];
                    elseif coord(bound.lnods(ib,1),2)-coord(bound.lnods(ib,2),2)>0;
                        aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                        aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    end
                    bound.normals = [bound.normals ; nleft'];
                    bound.normals = [bound.normals ; nleft'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                end
            end
        % The check below are for elements located at corner(s) and have 
        % nodes on 2 boundaries or more than 2 corners.
        % Take the bottom-right corner as example, besides the normal 
        % situations such as 234 22344, have been considered, exerteme 
        % situations also covered, such as 12345 12344 22345 134 135 235;
        % similar considerations been applied to all four corners
        elseif (~isequal((number == [16100 17200 26200 37100 46100 25100 36000 45000]),zeros(1,8)))
            if nnode == 3 || nnode == 4;
                for j=1:3;
                    if bound.nodpos(ib,j) == 1
                        j1=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 2
                        j1=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 3
                        j3=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 4
                        j5=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 5
                        j5=bound.lnods(ib,j);
                    end
                end
                % In this part, since the bottom right corner is fixed, as
                % long as the node is happened on the bottom boundary
                % (either on the boundary or the bottom left corner equally
                % bound.nodpos = 1 or 2) it will be considered as j1, same
                % situation for the nodes on the right boundary. In this
                % 'reduced' way aux will be more easily obtained
                aux1 = [aux1 ; j1 j3 ];
                bound.normals = [bound.normals ; nbottom'];
                aux1 = [aux1 ; j3 j5 ];
                bound.normals = [bound.normals ; nright'];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 6 || nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:5) , [2 2 3 4 4]);
                    % Again, since lnods might start from any corner, so
                    % doing such judgement is necessary.
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nright'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif isequal(bound.nodpos(ib,1:5) , [3 4 4 2 2]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif isequal(bound.nodpos(ib,1:5) , [4 2 2 3 4]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                end
            end
        % The top right corner
        elseif (~isequal((number == [20110 20220 70210 35120 85110 85000 70100 35010]),zeros(1,8)))
            if nnode == 3 || nnode == 4;
                for j=1:3;
                    if bound.nodpos(ib,j) == 3
                        j3=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 4
                        j3=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 5
                        j5=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 6
                        j7=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 7
                        j7=bound.lnods(ib,j);
                    end
                end
                aux1 = [aux1 ; j3 j5 ];
                bound.normals = [bound.normals ; nright'];
                aux1 = [aux1 ; j5 j7 ];
                bound.normals = [bound.normals ; ntop'];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 6 || nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:5) , [4 4 5 6 6]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; ntop'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif isequal(bound.nodpos(ib,1:5) , [5 6 6 4 4]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif isequal(bound.nodpos(ib,1:5) , [6 4 4 5 6]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                end
            end
        % The top left corner
        elseif (~isequal((number == [50011 50022 70012 60021 80011 60010 70001 80000]),zeros(1,8)))
            if nnode == 3 || nnode == 4;
                for j=1:3;
                    if bound.nodpos(ib,j) == 1
                        j1=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 5
                        j5=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 6
                        j5=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 7
                        j7=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 8
                        j1=bound.lnods(ib,j);
                    end
                end
                aux1 = [aux1 ; j5 j7 ];
                bound.normals = [bound.normals ; ntop'];
                aux1 = [aux1 ; j7 j1 ];
                bound.normals = [bound.normals ; nleft'];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 6 || nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:5) , [6 6 7 8 8]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nleft'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif isequal(bound.nodpos(ib,1:5) , [7 8 8 6 6]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif isequal(bound.nodpos(ib,1:5) , [8 6 6 7 8]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                end
            end
        % The bottom left corner
        elseif (~isequal((number == [11001 12002 26002 62001 76001 25001 61000 75000]),zeros(1,8)))
            if nnode == 3 || nnode == 4;
                for j=1:3;
                    if bound.nodpos(ib,j) == 1
                        j1=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 2
                        j3=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 3
                        j3=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 7
                        j7=bound.lnods(ib,j);
                    elseif bound.nodpos(ib,j) == 8
                        j7=bound.lnods(ib,j);
                    end
                end
                aux1 = [aux1 ; j7 j1 ];
                bound.normals = [bound.normals ; nleft'];
                aux1 = [aux1 ; j1 j3 ];
                bound.normals = [bound.normals ; nbottom'];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 6 || nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:5) , [1 2 2 8 8]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif isequal(bound.nodpos(ib,1:5) , [2 8 8 1 2]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                elseif isequal(bound.nodpos(ib,1:5) , [8 8 1 2 2]);
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                    aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                end
            end
            
        % Line 681 to line 1173 are some wired 'possible' check-points for  
        % 4 noded and 8 or 9 noded element situarions.
        elseif (~isequal((number == [61010 62021]),zeros(1,2)))
            if nnode == 4;
                if isequal(bound.nodpos(ib,1:4) , [1 2 6 7]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [2 6 7 1]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [6 7 1 2]);
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [7 1 2 6]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:7) , [1 2 2 6 6 7 8]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:7) , [2 6 6 7 8 1 2]);
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:7) , [6 6 7 8 1 2 2]);
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:7) , [7 8 1 2 2 6 6]);
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            end
        elseif (~isequal((number == [36010 37120]),zeros(1,2)));
            if nnode == 4;
                if isequal(bound.nodpos(ib,1:4) , [2 3 5 6]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:4) , [3 5 6 2]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                elseif isequal(bound.nodpos(ib,1:4) , [5 6 2 3]);
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:4) , [6 2 3 5]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:7) , [2 2 3 4 5 6 6]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:7) , [3 4 5 6 6 2 2]);
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:7) , [5 6 6 2 2 3 4]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:7) , [6 2 2 3 4 5 6]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            end
        elseif (~isequal((number == [25101 26202]),zeros(1,2)));
            if nnode == 4;
                if isequal(bound.nodpos(ib,1:4) , [1 3 4 8]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [3 4 8 1]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [4 8 1 3]);
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [8 1 3 4]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:7) , [1 2 3 4 4 8 8]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:7) , [3 4 4 8 8 1 2]);
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:7) , [4 8 8 1 2 3 4]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:7) , [8 8 1 2 3 4 4]);
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            end
        elseif (~isequal((number == [70110 70212]),zeros(1,2)));
            if nnode == 4;
                if isequal(bound.nodpos(ib,1:4) , [8 4 5 7]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [4 5 7 8]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [5 7 8 4]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [7 8 4 5]);
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:7) , [8 4 4 5 6 7 8]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:7) , [4 4 5 6 7 8 8]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:7) , [5 6 7 8 8 4 4]);
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:7) , [7 8 8 4 4 5 6]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,7)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,7) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            end
        elseif (~isequal((number == [2020 3030]),zeros(1,2)));
            if nnode == 4;
                if isequal(bound.nodpos(ib,1:4) , [2 2 6 6]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:4) , [2 6 6 2]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:4) , [6 6 2 2]);
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:4) , [6 2 2 6]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:6) , [2 2 2 6 6 6]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:6) , [2 6 6 6 2 2]);
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:6) , [6 6 6 2 2 2]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; ntop'];
                elseif isequal(bound.nodpos(ib,1:6) , [6 2 2 2 6 6]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nbottom'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; ntop'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; ntop'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            end
        elseif (~isequal((number == [202 303]),zeros(1,2)));
            if nnode == 4;
                if isequal(bound.nodpos(ib,1:4) , [8 4 4 8]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [4 4 8 8]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [4 8 8 4]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:4) , [8 8 4 4]);
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            elseif nnode == 8 || nnode == 9;
                if isequal(bound.nodpos(ib,1:6) , [8 4 4 4 8 8]);
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:6) , [4 4 4 8 8 8]);
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:6) , [4 8 8 8 4 4]);
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,6) bound.lnods(ib,1)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,3) bound.lnods(ib,4)];
                    bound.normals = [bound.normals; nleft'];
                elseif isequal(bound.nodpos(ib,1:6) , [8 8 8 4 4 4]);
                    aux1 = [aux1; bound.lnods(ib,4) bound.lnods(ib,5)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,5) bound.lnods(ib,6)];
                    bound.normals = [bound.normals; nright'];
                    aux1 = [aux1; bound.lnods(ib,1) bound.lnods(ib,2)];
                    bound.normals = [bound.normals; nleft'];
                    aux1 = [aux1; bound.lnods(ib,2) bound.lnods(ib,3)];
                    bound.normals = [bound.normals; nleft'];
                end
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
                aux3 = [aux3 ; bound.elno(ib)]; aux3 = [aux3 ; bound.elno(ib)];
            end
        end
    end
    % boundary face normals
    bound.normals = bound.normals(2:length(bound.normals),:);
    % boundary face connectivities
    bound.lnods   = aux1(2:length(aux1),:);
    % shows the element number of each boundary face
    bound.elno    = aux3(2:length(aux3));
    % Assemble integral constraint matrix
    Cglob = zeros(3,ndofn*npoin);
    
    % below are the solution for T3 type elements traction boundary constrain
    % condition
    if nnode == 3 || nnode == 4;
        for il = 1:length(bound.elno);
            % face length
            c1 = coord(bound.lnods(il,1),:); c2 = coord(bound.lnods(il,2),:);
            facelen = norm(c1-c2);
            % components of unit normal to face
            nx = bound.normals(il,1); ny = bound.normals(il,2);
            % constraint matrix (element-by-element)
            a = facelen*0.5;
            Celem(:,:,il) = a*[  nx       0      nx        0   ;
                0      ny       0       ny   ;
                0.5*ny  0.5*nx  0.5*ny  0.5*nx ];
            dof1 = bound.lnods(il,1)*ndofn-1;
            dof2 = dof1+1;
            dof3 = bound.lnods(il,2)*ndofn-1;
            dof4 = dof3+1;
            dofglob = [ dof1 dof2 dof3 dof4 ];
            %...global constraint matrix
            Cglob(:,dofglob) = Cglob(:,dofglob)+Celem(:,:,il);
        end
        Cd(:,:) = Cglob(:,dofs.d);
        Cf(:,:) = Cglob(:,dofs.f);
        bound.R = -inv(Cd)*Cf;
    elseif nnode == 6 || nnode == 8 || nnode == 9;
        for il = 1:length(bound.elno)*0.5;
            % face length
            c1 = coord(bound.lnods(2*il-1,1),:); c2 = coord(bound.lnods(2*il,2),:);
            facelen = norm (c1-c2);
            % components of unit normal to face
            nx = bound.normals(2*il-1,1); ny = bound.normals(2*il-1,2);
            % use the following check function to determin the element
            % boundary on the horizontal side (top & bottom) or on the
            % vertical side (left & right), so that could assign numbers to
            % x1, x2, x3, and carry out gauss integration
            if coord(bound.lnods(2*il-1,1),1)-coord(bound.lnods(2*il,2),1)==0;
                x1=coord(bound.lnods(2*il-1,1),2); x2=coord(bound.lnods(2*il-1,2),2);
                x3=coord(bound.lnods(2*il,2),2);
            else
                x1=coord(bound.lnods(2*il-1,1),1); x2=coord(bound.lnods(2*il-1,2),1);
                x3=coord(bound.lnods(2*il,2),1);
            end
            if x1 > x3;
                x4 = x1;
                x1 = x3;
                x3 = x4;
            end
            % due to the heighest order of the shape function is 2, so 2
            % gauss points is accurate enough to do the calculation due to
            % the weights for 2 gauss points are 1, so just negelect
            posgp(1) = -5.773502691896258e-001; % -1/sqrt(3)
            posgp(2) = 5.773502691896258e-001;
            % define the new 'x' in gauss quadrature as z1 and z2
            z1 = (x1+x3)/2+(x3-x1)/2*posgp(1);
            z2 = (x1+x3)/2+(x3-x1)/2*posgp(2);
            a1 = (facelen/2)*(((z1-x2)*(z1-x3))+((z2-x2)*(z2-x3)))/((x1-x2)*(x1-x3));
            a2 = (facelen/2)*(((z1-x1)*(z1-x3))+((z2-x1)*(z2-x3)))/((x2-x1)*(x2-x3));
            a3 = (facelen/2)*(((z1-x1)*(z1-x2))+((z2-x1)*(z2-x2)))/((x3-x2)*(x3-x1));
            % constraint matrix (element-by-element)
            Celem(:,:,il) = [  a1*nx     0     a2*nx     0     a3*nx     0;
                0     a1*ny     0     a2*ny     0     a3*ny;
                0.5*a1*ny 0.5*a1*nx  0.5*a2*ny  0.5*a2*nx 0.5*a3*ny 0.5*a3*nx];
            dof1 = bound.lnods(2*il-1,1)*ndofn-1;
            dof2 = dof1+1;
            dof3 = bound.lnods(2*il-1,2)*ndofn-1;
            dof4 = dof3+1;
            dof5 = bound.lnods(2*il,2)*ndofn-1;
            dof6 = dof5+1;
            dofglob = [ dof1 dof2 dof3 dof4 dof5 dof6];
            %...global constraint matrix
            Cglob(:,dofglob) = Cglob(:,dofglob)+Celem(:,:,il);
        end
        Cd(:,:) = Cglob(:,dofs.d);
        Cf(:,:) = Cglob(:,dofs.f);
        bound.R = -inv(Cd)*Cf;
    end
end
end