function [ke] = calcLocalSystem(eldof,NGP,D,nen,eltype,x0,elc,ie,GQ)
       
% initialize the submatrices and vectors of the local system
ke = zeros(eldof,eldof);
B = zeros(3,eldof);  % for 3D problem B = zeros(6,eldof)
                   

% Generate nodal coordinate matrix of size obj.nen x 2. Each row of it 
% stores x and y coords of the nodes of elem e. It is used to calculate 
% the Jacobian matrix.
e_x0 = x0(elc(:,ie),:);
                

% Calculate elemental matrix.
for k = 1:NGP    % Gauss points loop 
    
    % Calculate shape functions and Jacobian matrix.
    [~,gdN,Jacob] = calcShapefunAndJacob(k,nen,eltype,GQ.point,e_x0);
    
    % Define shortcuts.
    factor = det(Jacob) * GQ.weight(k);    
    
    
    for i = 1:nen
        i1 = (i-1) * 2 + 1;
        i2 = i1 + 1;
        B(1,i1) = gdN(1,i);
        B(2,i2) = gdN(2,i);
        B(3,i1) = gdN(2,i);
        B(3,i2) = gdN(1,i);
    end
    ke = ke + B' * D * B * factor;
    
end  % End of GQ loop