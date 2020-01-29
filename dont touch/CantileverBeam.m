clc;clear all;close all;
%
% Set some basic variables
% ------------------------
% Number of spatial dimensions
ndim = 2;
% Number of degrees of freedom per node
ndofn = 2;
%
%==================================================
%   DATA INPUT PHASE. READ DATA FROM DATA FILE
%==================================================
%
fid         = fopen('CantileverQua.inp');
title       = fscanf(fid,'TITLE = %s',1);
disp('');disp('');
disp('==================================================================');
disp(['Problem title: ',title]);
%
% Define number of element groups
%
nelgroups   = fscanf(fid,'\nNUMBER OF ELEMENT GROUP = %d',1);
groupcheck  = zeros(nelgroups,1);
%
% READ MATERIAL PROPERTIES FOR ALL ELEMENT GROUPS
%
for ielgroup = 1:nelgroups
    groupid = fscanf(fid, '\nELEMENT_GROUP = %d', 1);
    if groupid <=0 | groupid > nelgroups
        error('ELEMENT_GROUP number smaller or equal to zero or greater than specified number of material groups')
    end
    if groupcheck(groupid)
        error('Repeated ELEMENT_GROUP number specified in data file')
    else
        groupcheck(groupid)=1;
    end
    matprop(groupid).young = fscanf(fid,'\nYOUNG_MODULUS = %E',1);
    matprop(groupid).poiss = fscanf(fid,'\nPOISSON_RATIO = %f',1);
end
%
% Read element properties according to specified element type
%
elemtype    = fscanf(fid, '\nELEMENT_TYPE = %s',1);
if strcmp(elemtype,'TRI_3')
    % Three-Node (simplex) triangle
    nnode = 3;
    ngaus = fscanf(fid,'\nNUMBER_OF_GP = %d',1);
    domain = 'TRI';
elseif strcmp(elemtype,'TRI_6')
    % six-Noded Quadratic quadrilateral
    nnode = 6;
    ngaus = fscanf(fid,'\nNUMBER_OF_GP = %d',1);
    domain = 'TRI';
elseif strcmp(elemtype,'QUAD_4')
    % four-Noded Quadratic quadrilateral
    nnode = 4;
    ngaus = fscanf(fid,'\nNUMBER_OF_GP = %d',1);
    domain = 'QUA';
elseif strcmp(elemtype,'QUAD_8')
    % Eight-noded quadratic quadrilateral
    nnode = 8;
    ngaus=fscanf(fid,'\nNUMBER_OF_GP = %d', 1);
    domain = 'QUA';
elseif strcmp(elemtype,'QUAD_9')
    % Nine-noded quadratic quadrilateral
    nnode = 9;
    ngaus=fscanf(fid,'\nNUMBER_OF_GP = %d', 1);
    domain = 'QUA';
end
%
%   Read analysis type
%
%   ntype = 1 -> Plane Stress Analysis
%   ntype = 2 -> Plane Strain analysis
%
ntype = fscanf (fid,'\nANALYSIS_TYPE = %d',1);
disp('-------------------------------------------------------------------')
if ntype ==1
    disp ('Plane stress assumption');
elseif ntype ==2
    disp ('Plane strain assumption');
else
    error('Invalid ANALYSIS_TYPE specification in data file.Valid values are 1 (plane stress) and 2 (Plane strain).');
end
disp('-------------------------------------------------------------------')
%
% Read RVE mesh data
%
%...total number of elements
nelem = fscanf (fid,'\nELEMENTS = %d',1);
%
%...table of connectivities and element group information
%
lnods = fscanf(fid,'\n%d', [2+nnode, nelem]);
lnods = lnods';
elorder = lnods(:,1);
lnods(elorder,:) = lnods;
%
% build up the lnods_final
%
lnodes_final = zeros(nelem,nnode);
elgroup = lnods(:,2);
if sum(elgroup<=0)~=0 || sum(elgroup>nelgroups)~=0
    error('At least one element has been assigned an invalid group number')
end
lnods = lnods(:,2:nnode+1);
%...store the set of element numbers of each element group separately
%   (used for plotting purposes only)
groupel = zeros(1,nelgroups);
for ielgroup = 1:nelgroups
    count = 0;
    for ielem = 1:nelem
        if elgroup(ielem) == ielgroup
            count = count +1;
            groupel(count,ielgroup) = ielem;
        end
    end
    nelemg(ielgroup) = count;
end
%
%...nodal coordinates
%
npoin = fscanf(fid, '\nNODE_COORDINATES = %d', 1);
coord = fscanf(fid, '\n%d %d %d', [3, npoin]);
coord = coord';
coord = coord(:,2:3);
%
%...material thickness (plane stress only)
%
if ntype == 1
    thick = fscanf(fid, '\nTHICKNESS = %f', 1);
else
    thick = 0;
end
% Close data file before proceeding to solution phase
%
status = fclose(fid);
% =====================================================================
%                          SOLUTION PHASE
% =====================================================================
%
% Compute global stiffness matrix
%
K = zeros(npoin*ndofn,npoin*ndofn);
if ndim == 2
    B = zeros(3,npoin*ndofn,nelgroups);
end
for ielgroup = 1:nelgroups
    % Retrieve material properties
    young = matprop.young; poiss = matprop.poiss;
    %
    % Compute elasticity tensor in matrix form (D matrix)
    %
    if  ntype == 1      % plane stress
        const(ielgroup) = young/(1-poiss^2);
        dmatx(:,:,ielgroup) = const(ielgroup)*[  1    poiss       0       ;
                                               poiss    1         0       ;
                                                 0      0    (1-poiss)/2 ];
    elseif  ntype == 2  % Plane strain
        const(ielgroup) = young/((1+poiss)*(1-2*poiss));
        dmatx(:,:,ielgroup) = const(ielgroup)*[1-poiss   poiss          0        ;
                                                poiss   1-poiss         0        ;
                                                   0        0       (1-2*poiss)/2];
    end
    for ielem = groupel(:,ielgroup)'
        if ielem == 0
            break
        end
        % Compute element stiffness and B matrix
        [area(ielem),ke(:,:,ielem),intbmatx(:,:,ielem),lnods_trans(:,:,ielem),accbmatx(:,:,ielem)] = stiffness(dmatx,domain,ngaus,coord,ielem,lnods,ntype,thick,nnode,ndim,ielgroup);
        lnods_final(ielem,:) = lnods_trans(:,:,ielem);
        % assemble table of global dofs of current element
        clear gdof;
        ildof = 0;
        for inode = 1:nnode
            for idofn = 1:ndofn
                ildof = ildof + 1;
                gdof(ildof) = (lnods_final(ielem,inode)-1)*ndofn+idofn;
            end
        end
        % add contribution of current element stiffness to global
        % stiffness
        K(gdof(:),gdof(:)) = K(gdof(:),gdof(:)) + ke(:,:,ielem);
        % add contribution of current element to the product B(global)*vol
        B(:,gdof(:),ielgroup) = B(:,gdof(:),ielgroup) + intbmatx(:,:,ielem);
    end
end
%
%
% Computation of Homogenised Elasticity Tensor
% --------------------------------------------
%
% Compute G Global Matrix
if ndim == 2
    G = zeros(3,npoin*ndofn);
end
for ielgroup = 1:nelgroups
    young = matprop(ielgroup).young; poiss = matprop(ielgroup).poiss;
    %
    % Compute elasticity tensor in matrix form (D matrix)
    %
    if  ntype == 1      % plane stress
        const = young/(1-poiss^2);
        dmatx(:,:,ielgroup) = const*[  1    poiss       0       ;
                                     poiss    1         0       ;
                                       0      0    (1-poiss)/2 ];
    elseif  ntype == 2  % Plane strain
        const = young/((1+poiss)*(1-2*poiss));
        dmatx(:,:,ielgroup) = const*[1-poiss   poiss          0        ;
                                      poiss   1-poiss         0        ;
                                        0        0       (1-2*poiss)/2];
    end
    G(:,:) = G(:,:) + dmatx(:,:,ielgroup)*B(:,:,ielgroup);
end
G=G';
% Split mesh nodes and dof's into convenient sets
[bound,cell_volume,dofs,nodes] = split_mesh(coord,lnods_final,nelem,ndofn,npoin,option,nnode);
%
% Compute Taylor Constitutive Matrix
%
Dtaylor=[[0 0 0];
    [0 0 0];
    [0 0 0]];
volfractot=0;
for ielgroup = 1:nelgroups
    matpropg = matprop(ielgroup);
    %... Micro-Elasticity Matrix for current phase (group)
    Du=dmatx(:,:,ielgroup);
    indices=groupel(:,ielgroup)~=0;
    vol=sum(area(groupel(indices,ielgroup)));
    volfrac = vol/cell_volume;volfractot=volfractot+volfrac;
    Dtaylor = Dtaylor + volfrac*Du;
    disp(['Phase no. ',num2str(ielgroup),':   E = ',num2str(matprop(ielgroup).young),'  nu = ',num2str(matprop(ielgroup).poiss),'  Vol.frac. = ',num2str(volfrac)]);
end
disp(['Total vol.frac. of solid phases: ',num2str(volfractot)]);
disp('------------------------------------------------------------------')
% Compute Homogenised Constitutive Matrix
if strcmp(option,'taylor')
    disp('RVE kinematics: Taylor')
    Dhom=Dtaylor;
elseif strcmp(option,'linear')
    disp('RVE kinematics: Linear boundary displacements')
    disp('------------------------------------------------------------------')
    Kii = K(dofs.i(:),dofs.i(:));
    Kr=Kii;
    Gi = G(dofs.i(:),:);
    Dhom = Dtaylor - (1/cell_volume)*Gi'/Kr*Gi;
elseif strcmp(option,'periodic')
    disp('RVE kinematics: Periodic boundary displacement fluctuations')
    disp('------------------------------------------------------------------')
    % Set relevant stiffness sub-matrices
    Kii = K(dofs.i(:),dofs.i(:));
    Kip = K(dofs.i(:),dofs.p(:));
    Kim = K(dofs.i(:),dofs.m(:));
    Kpi = K(dofs.p(:),dofs.i(:));
    Kpp = K(dofs.p(:),dofs.p(:));
    Kpm = K(dofs.p(:),dofs.m(:));
    Kmi = K(dofs.m(:),dofs.i(:));
    Kmp = K(dofs.m(:),dofs.p(:));
    Kmm = K(dofs.m(:),dofs.m(:));
    % Assemble reduced system matrix
    Kr = [   Kii         Kip+Kim      ;
        Kpi+Kmi   Kpp+Kpm+Kmp+Kmm ];
    Gi = G(dofs.i(:),:);
    Gp = G(dofs.p(:),:);
    Gm = G(dofs.m(:),:);
    Gr = [    Gi    ;
        Gp+Gm  ];
    Dhom = Dtaylor - (1/cell_volume)*Gr'/Kr*Gr;
elseif strcmp(option,'traction')
    disp('RVE kinematics: Minumum RVE constraint (uniform boundary traction)')
    disp('------------------------------------------------------------------')
    % Set relevant stiffness sub-matrices
    Kii = K(dofs.i(:),dofs.i(:));
    Kif = K(dofs.i(:),dofs.f(:));
    Kid = K(dofs.i(:),dofs.d(:));
    Kfi = K(dofs.f(:),dofs.i(:));
    Kff = K(dofs.f(:),dofs.f(:));
    Kfd = K(dofs.f(:),dofs.d(:));
    Kdi = K(dofs.d(:),dofs.i(:));
    Kdf = K(dofs.d(:),dofs.f(:));
    Kdd = K(dofs.d(:),dofs.d(:));
    % Retrieve integral constraint matrix
    R = bound.R;
    % Assemble reduced system stiffness
    Kr = [   Kii               Kif+Kid*R           ;
        Kfi+R'*Kdi  Kff+Kfd*R+R'*Kdf+R'*Kdd*R  ];
    Gi = G(dofs.i(:),:);
    Gf = G(dofs.f(:),:);
    Gd = G(dofs.d(:),:);
    Gr = [     Gi      ;
        Gf+R'*Gd  ];
    Dhom = Dtaylor - (1/cell_volume)*Gr'/Kr*Gr;
end
disp(' ');
disp('Homogenised elasticity tensor (in matrix format):');disp(' ');
disp(num2str(Dhom));disp(' ');
%
% Compute the homogenised stress tensor
% -------------------------------------
if stressout
    % Engineering strain
    estrain=[strain(1); strain(2); 2*strain(3)];
    % Homogenised stress
    stress=Dhom*estrain;
    disp('------------------------------------------------------------------')
    disp(' ');
    disp('Homogenised strains and stresses:');disp(' ');
    disp([['exx = ';'eyy = ';'exy = '],num2str(estrain),['    sxx = ';'    syy = ';'    sxy = '],num2str(stress)]);disp(' ');
end
disp('==================================================================')
disp(' ');disp(' ');
%
% Compute RVE displacement and stress fields
% ------------------------------------------
if stressout
    strain_tensor = [  strain(1)    strain(3)   ;
                       strain(3)    strain(2)  ];
    % Set prescribed component (uniform strain component) of all
    % nodal displacements of RVE mesh
    % ...i nodes
    Ustar = [];
    for i = 1:length(nodes.all)
        y  = coord(nodes.all(i),:)';
        Ustar = [Ustar; strain_tensor*y ];
    end
    if strcmp(option,'taylor')
        % All displacement fluctuations are prescribed as zero
        Utilde = zeros(length(nodes.all)*2,1);
    elseif strcmp(option,'linear')
        % Assemble reduced system forcing term
        Kstar = K(dofs.i(:),dofs.all(:));
        Fr = - Kstar * Ustar;
        % Solve reduced system for unknwon displacement fluctuations
        % (interior nodes only)
        Ur = Kr \ Fr;
        Uitilde = Ur;
        %... all boundary node fluctuations fully constrained
        Ubtilde = zeros(length(nodes.b)*2,1);
        %... assemble vector of nodal displ fluctuations
        Utilde = [ Uitilde ; Ubtilde ];
    elseif strcmp(option,'periodic')
        % Set other relevant stiffness sub-matrices
        Kic = K(dofs.i(:),dofs.c(:));
        Kpc = K(dofs.p(:),dofs.c(:));
        Kmc = K(dofs.m(:),dofs.c(:));
        % Assemble reduced system forcing term
        Kstar = [    Kii      Kip      Kim      Kic     ;
                   Kpi+Kmi  Kpp+Kmp  Kpm+Kmm  Kpc+Kmc  ];
        Fr = - Kstar * Ustar;
        % Solve reduced system for unknown periodic displacements
        Ur = Kr \ Fr;
        ngdof_interior = length(dofs.i);
        Uitilde = Ur(1:ngdof_interior);
        Uptilde = Ur(ngdof_interior+1:length(Ur));
        Umtilde = Uptilde;
        %... all corner nodes fully constrained
        Uctilde = zeros(length(nodes.c)*2,1);
        %... assemble vector of nodal displ fluctuations
        Utilde = [ Uitilde ; Uptilde ; Umtilde ; Uctilde ];
    elseif strcmp(option,'traction')
        % Set other relevant stiffness sub-matrices
        Kiz = K(dofs.i(:),dofs.z(:));
        Kfz = K(dofs.f(:),dofs.z(:));
        Kdz = K(dofs.d(:),dofs.z(:));
        % Assemble reduced system forcing term
        Kstar = [     Kii         Kif         Kid         Kiz      ;
                   Kfi+R'*Kdi  Kff+R'*Kdf  Kfd+R'*Kdd  Kfz+R'*Kdz ];
        Fr = - Kstar * Ustar;
        % Solve reduced system for independent dof's
        Ur = Kr \ Fr;
        ni = length(dofs.i); nf = length(dofs.f);
        Uitilde = Ur(1:ni);
        Uftilde = Ur(ni+1:ni+nf);
        %... compute dependent dof's
        Udtilde = R*Uftilde;
        %... assemble vector of nodal displ fluctuations
        Utilde = [ Uitilde ; Uftilde ; Udtilde ; zeros(3,1) ];
    end
    %
    % Compute full nodal displacement vector
    %
    U = Ustar + Utilde;
    % Assemble micro-cell global nodal displacement vector
    % (use global nodal numbering)
    Ug([dofs.all]) = U;
    %...split U into x and y components for plotting
    Ux = Ug(1:2:length(Ug))'; Uy = Ug(2:2:length(Ug))';
end
%
% =====================================================================
% POST-PROCESSING PHASE. PLOT RVE MESH (AND EFFECTIVE STRESS FIELD IF
% STRAINS ARE GIVEN IN INPUT ARGUMENT LIST)
% =====================================================================
%
hold off;clf; xlim('auto');ylim('auto');
%
% Plot RVE mesh (undeformed state)
% --------------------------------
%
% Split mesh into groups to plot mesh using a different color for each
% element group
colorarray = [ 'r' ; 'b' ; 'g' ; 'm' ; 'k' ; 'c' ; 'y' ];
ncolor = length(colorarray);
for ielgroup = 1:nelgroups
    eleconvert = groupel(:,ielgroup);
    nodno  = zeros(npoin,1);
    for ielemg = 1:nelemg(ielgroup)
        iel = groupel(ielemg,ielgroup);
        eleconvert(ielemg) = iel;
        for inode = 1:nnode
            nodno(lnods(iel,inode)) = 1;
        end
    end
    %
    npoing = 0;
    nodconvert  = zeros(npoin,1);
    nodconvert2 = zeros(npoin,1);
    for ipoin = 1:npoin
        if nodno(ipoin) == 1
            npoing = npoing + 1;
            nodconvert(npoing) = ipoin;
            nodconvert2(ipoin) = npoing;
        end
    end
    if nnode == 9;
        coordg = coord(nodconvert(1:npoing),1:8);
    else
        coordg = coord(nodconvert(1:npoing),:);
    end
    % assemble connectivity table and plot current element group
    lnodsg = nodconvert2(lnods(eleconvert(1:nelemg(ielgroup)),1:nnode));
    if ielgroup <= ncolor
        color = colorarray(ielgroup,:);
    else
        color = '-k';
    end
    patch('Faces',lnodsg,'Vertices',coordg,'FaceColor','w','EdgeColor',color);
    hold on;
end
%
% Plot deformed mesh with von Mises effective stress contour (only if a
% "strain" argument is passed into this function)
% ---------------------------------------------------------------------
%
if stressout
    % Calculate the total area that a node assigned
    for ip = 1:npoin;
        area_total = 0;
        for ielem = 1:nelem;
            if (~isequal(ip == lnods_final(ielem,:),zeros(1,nnode)));
                area_total = area_total + area(ielem);
            end
            area_final(ip) = area_total;
        end
    end
    nodeffstr=zeros(npoin,1);nodestrs=zeros(npoin,3);Uel=zeros(nelem,nnode*2);
    for ielgroup = 1:nelgroups
        % Compute stress in each element
        for ielem = groupel(:,ielgroup)'
            if ielem == 0
                break
            end
            %...assemble table of global dofs of current element
            clear gdof;
            ildof = 0;
            for inode = 1:nnode
                for idofn = 1:ndofn
                    ildof = ildof + 1;
                    gdof(ildof) = (lnods_final(ielem,inode)-1)*ndofn+idofn;
                end
            end
            %...retrieve element displacements
            Uel(ielem,:)=Ug(gdof);
            [posgp] = gaus2d(domain,ngaus);
            for igaus = 1:ngaus;
                s=posgp(1,igaus);
                t=posgp(2,igaus);
                [deriv] = SFandD(s,t,nnode);
                for in = 1:nnode;
                    eloc(in,:)=coord(lnods_final(ielem,in),:);
                end
                jac = deriv*eloc;
                cartd = inv(jac)*deriv;
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
                % Calculate the stress at each gauss point in each
                % element from each element group
                strsg(:,igaus,ielem) = dmatx(:,:,ielgroup)*bmatx*Uel(ielem,:)';
            end
            [EXMATX] = EXPLT_F(nnode,ngaus);
            elastra(:,:,ielem) = EXMATX*strsg(:,:,ielem)';
            for io = 1:nnode;
                nodestrs(lnods_final(ielem,io),:) = nodestrs(lnods_final(ielem,io),:)+elastra(io,:,ielem)*area(ielem);
            end
        end
    end
    %...add contribution of current element to its nodal values of
    %   smoothed effective stress
    %...area averaged smoothed effective stress
    for iq = 1:npoin;
        nodestrs(iq,:) = nodestrs(iq,:)./area_final(iq)';
    end
    clear ielem
    clear for ielgroup = 1:nelgroups
    for ielgroup = 1:nelgroups
        for ielem = 1:nelem;
            for is = 1:nnode;
                if ntype==1
                    % plane stress
                    elstres_tensor=[ nodestrs(lnods_final(ielem,is),1)  nodestrs(lnods_final(ielem,is),3) 0  ;
                        nodestrs(lnods_final(ielem,is),3)  nodestrs(lnods_final(ielem,is),2) 0  ;
                        0            0      0 ];
                elseif ntype==2
                    % plane strain
                    elstra = accbmatx(:,:,ielem)*Uel(ielem,:)';
                    szz=const*matprop(elgroup(ielem)).poiss*(elstra(1)+elstra(2));
                    elstres_tensor=[ nodestrs(lnods_final(ielem,is),1)  nodestrs(lnods_final(ielem,is),3)   0  ;
                        nodestrs(lnods_final(ielem,is),3)  nodestrs(lnods_final(ielem,is),2)   0  ;
                        0           0        szz];
                end
                devstres=elstres_tensor-trace(elstres_tensor)/3*[1 0 0;0 1 0;0 0 1];
                %...von Mises effective stress (smoothed - area averaged)
                nodeffstr(lnods_final(ielem,is))=sqrt(3/2*sum(sum(devstres.*devstres)));
            end
        end
    end
    %
    % Plot effective stress contours on deformed mesh
    %...apply rigid displacement to initial mesh
    rdisp = [max(coord(:,1))-min(coord(:,1)) max(coord(:,2))-min(coord(:,2))]*1.1;
    %...compute  deformed coordinates (with added rigid displacements, so
    %   as to not overlap with underformed mesh plot, and amplified
    %   displacements to ease visualisation)
    defcoord(:,1) = coord(:,1) + amplfact*Ux + rdisp(1);
    defcoord(:,2) = coord(:,2) + amplfact*Uy + rdisp(2);
    if nnode == 9;
        lnods_final = lnods_final(:,1:8);
    end
    patch('Faces',lnods_final,'Vertices',defcoord,'FaceVertexCData',nodeffstr,'FaceColor','interp','EdgeColor','k');
    colorbar('East');
    
end
%
% Re-set plot limits
xmax=max(xlim);xmin=min(xlim);ymax=max(ylim);ymin=min(ylim);
lx = xmax-xmin; ly = ymax-ymin;
lplot = max([lx;ly]);
xlim([xmin-lplot*0.18 xmin+lplot*1.18]);ylim([ymin-lplot*0.05 ymin+lplot*1.05])
