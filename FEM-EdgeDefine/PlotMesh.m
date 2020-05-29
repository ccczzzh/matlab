function PlotMesh(nodecoord,node,ne,nn,nen,show)
if nargin == 5
    show = 0;
end
dimension = size(nodecoord,2);

% 
% Initialization of the required matrices
X = zeros(nen,ne) ;
Y = zeros(nen,ne) ;
% Z = zeros(nen,ne) ;

if dimension == 2
    for i = 1:ne
        nd = node(i,:);
        X(:,i) = nodecoord(nd,1);
        Y(:,i) = nodecoord(nd,2);
    end
    
    figure; fill(X,Y,'w');
    title('Finite Element Mesh')
    axis off;
    if show ~= 0
        k = 1:nn;
        nd = k';
        for i = 1:ne
            text(X(:,i),Y(:,i),int2str(nd(i)),'fontsize',8,'color','k');
            text(sum(X(:,i))/4,sum(Y(:,i))/4,int2str(i),'fontsize',10,'color','r') ;
        end
    end
    
end