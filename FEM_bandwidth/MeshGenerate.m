clc;clear;close;

length = input('Enter the length of the component:');
%length = 20;
height = input('Enter the height of the component:');
% height = 5;

seedx = input('Enter the number of the seed placed in the horizontal direction:');
% seedx = 5;
seedy = input('Enter the number of the seed placed in the vertical direction:');
% seedy = 3;

%eltype = input('Enter the type of elements: (3: triangular, 4: quadratical)');
eltype = 3;
Gcoord = zeros(seedx*seedy,2);
for ix = 1:(seedx-1)
    x_coord(ix+1) = ix*length/seedx;
    
end
for iy = 1:(seedy-1)
    y_coord(iy+1) = iy*height/seedy;
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
        nodemap((inen*2-1),:) = [inen, inen+seedx+1, inen+seedx];
        nodemap((inen)*2,:) = [inen, inen+1, inen+seedx+1];
        nodemap((inen*2-1+(ne/2)),:) = [inen+seedx, inen+seedx+seedx+1, inen+seedx+seedx];
        nodemap((inen)*2+(ne/2),:) = [inen+seedx, inen+seedx+1, inen+seedx+seedx+1];
    end
end
PlotMesh(Gcoord,nodemap,ne,nn,eltype)