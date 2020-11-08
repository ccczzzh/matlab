function PlotdeformMesh(Unew,node,ne,nen)
dimension = size(Unew,2);

X = zeros(nen,ne); Y=zeros(nen,ne);

if dimension == 2
    for i = 1:nen
        nd = node(i,:);
        X(i,:) = Unew(nd,1);
        Y(i,:) = Unew(nd,2);
    end
    
    figure(2);fill(X,Y,'w');
    title('deform mesh')
    axis off;
    
end