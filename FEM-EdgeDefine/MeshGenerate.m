function [nn,ne,nodemap,Gcoord] = MeshGenerate(seedx,seedy,eltype,ne,nn,Gcoord,LM)
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

% for iln = 1:localno
% Gcoord = zeros(seedx(iln)*seedy(iln),2);
% for ix = 1:(seedx(iln)-1)
%     x_coord(ix+1) = ix*length/(seedx(iln)-1);
%
% end
% for iy = 1:(seedy(iln)-1)
%     y_coord(iy+1) = iy*height/(seedy(iln)-1);
% end
%
% nnnew = seedx(iln)*seedy(iln);
% totalnn= 1:1:nnnew;
%
% % Node Coordinate
% for i = 1:seedy(iln)
%     for j = 1:seedx(iln)
%         n = j+seedx(iln)*(i-1);
%         jseedx = x_coord(j);iseedy = y_coord(i);
%        Gcoord(totalnn(n),:) = [jseedx, iseedy];
%     end
% end


if eltype == 3
    % Connectivity
    nenew = (seedx(iln)-1)*2*(seedy(iln)-1); % total number of elements - 16
    nodemap = zeros(nenew,eltype);% -[16,3]
    for inen = 1: (nenew/(2*(seedy(iln)-1)))
        for layer = 1:(seedy(iln)-1)
            nodemap((inen*2-1+(layer-1)*nenew/(seedy(iln)-1)),:) = [inen+(layer-1)*seedx(iln), inen+seedx(iln)+1+(layer-1)*seedx(iln), inen+seedx(iln)+(layer-1)*seedx(iln)];
            nodemap((inen*2+(layer-1)*nenew/(seedy(iln)-1)),:) = [inen+(layer-1)*seedx(iln), inen+1+(layer-1)*seedx(iln), inen+seedx(iln)+1+(layer-1)*seedx(iln)];
        end
    end
elseif eltype == 4
    %     ne = (seedx(iln)-1)*(seedy(iln)-1);
    % connectivity
    nodemap =zeros(ne,eltype);
    for inen = 1:(ne/(seedy-1))
        for layer = 1:(seedy-1)
            nodemap(inen+(layer-1)*(seedx-1),:) = [inen+(layer-1)*seedx, inen+1+(layer-1)*seedx, inen+seedx+1+(layer-1)*seedx,...
                inen+seedx+(layer-1)*seedx];
        end
    end
    %for iln = i:localno
    %        lrn = localrefineno(iln);
    %        localx = seedx(iln);localy = seedy(iln);
    %        [nodemap,localnodemap,nodemapreduce,Newnodemap,Lcoord,NewGcoord]=localrefinemesh(Gcoord,nodemap,nenew,nnnew,eltype,lrn,localx,localy,ndim);
    %    %end
end









