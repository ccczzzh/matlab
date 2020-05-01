function [nn,ne,nodemap,Gcoord] = MeshGenerate(length,height,seedx,seedy,eltype)
% length = input('Enter the length of the component:');
% 
% height = input('Enter the height of the component:');
% 
% seedx = input('Enter the number of the seed placed in the horizontal direction:');
% 
% seedy = input('Enter the number of the seed placed in the vertical direction:');

% seedy = 4;
% seedx = 5;
% height = 5;
% length = 20;

%eltype = input('Enter the type of elements: (3: triangular, 4: quadratical)');
% eltype = 4;
Gcoord = zeros(seedx*seedy,2);
for ix = 1:(seedx-1)
    x_coord(ix+1) = ix*length/(seedx-1);
    
end
for iy = 1:(seedy-1)
    y_coord(iy+1) = iy*height/(seedy-1);
end

nn = seedx*seedy;
totalnn= 1:1:nn;

% Node Coordinate
for i = 1:seedy
    for j = 1:seedx
        n = j+seedx*(i-1);
        jseedx = x_coord(j);iseedy = y_coord(i);
       Gcoord(totalnn(n),:) = [jseedx, iseedy];
    end
end


if eltype == 3
    % Connectivity
    ne = (seedx-1)*2*(seedy-1); % total number of elements - 16
    nodemap = zeros(ne,eltype);% -[16,3]
    for inen = 1: (ne/(2*(seedy-1)))
        for layer = 1:(seedy-1)
            nodemap((inen*2-1+(layer-1)*ne/(seedy-1)),:) = [inen+(layer-1)*seedx, inen+seedx+1+(layer-1)*seedx, inen+seedx+(layer-1)*seedx];
            nodemap((inen*2+(layer-1)*ne/(seedy-1)),:) = [inen+(layer-1)*seedx, inen+1+(layer-1)*seedx, inen+seedx+1+(layer-1)*seedx];
        end
    end
elseif eltype == 4
    ne = (seedx-1)*(seedy-1);
    % connectivity
   nodemap =zeros(ne,eltype);
   for inen = 1:(ne/(seedy-1))
      for layer = 1:(seedy-1)
         nodemap(inen+(layer-1)*(seedx-1),:) = [inen+(layer-1)*seedx, inen+1+(layer-1)*seedx, inen+seedx+1+(layer-1)*seedx, inen+seedx+(layer-1)*seedx]; 
      end
   end
end
PlotMesh(Gcoord,nodemap,ne,nn,eltype)









