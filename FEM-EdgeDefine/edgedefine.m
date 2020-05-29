%           edge 2
%         —————————
%        |              |
% edge 1 |              |   edge 3
%        |              |
%        |              |
%         --------------
%           edge 4
% clc;clear all;close all;
function [Gcoord,dn,node] = edgedefine(edge,X,Y,value,Gcoord,seedx,seedy)
seedx = 4 ; seedy = 5; 
if edge == 1
    node = zeros(seedy,1);
    for i = 1:seedy
        node(i) = 1 + (i-1)*seedx;
    end
elseif edge == 2
    node = zeros(seedx,1);
    for i = 1:seedx
        node(i) = i + seedx*(seedy-1);
    end
elseif edge == 3
    node = zeros(seedy,1);
    for i = 1:seedy
        node(i) = i*seedx;
    end
elseif edge == 4
    node = zeros(seedx,1);
    for i = 1:seedx
        node(i) = i;
    end
else
    error('Please Enter the edge number 1 - 4')
end

[n,m]=size(node);
% nd = n*2;nl = n;
nodecoord = zeros(n,2);
for i = 1:n
    nodecoord(i,:) = [node(i)*2-1,node(i)*2];
end
if X == 1
    Gcoord(nodecoord(:,1)) = value;
end
 if Y == 1
    Gcoord(nodecoord(:,2)) = value;
 end
dn(1:n,1) = nodecoord(:,1);
dn(n+1:2*n,1) = nodecoord(:,2);
dn = sort(dn);
