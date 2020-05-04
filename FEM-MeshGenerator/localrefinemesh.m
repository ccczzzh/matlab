function[nodemap,localnodemap,nodemapreduce,Newnodemap,Lcoord,NewGcoord] = localrefinemesh(Gcoord,nodemap,ne,nn,eltype,localrefineno,localx,localy,ndim)
% Local coordinate matrix.
Lcoord = zeros(localx*localy,2);
x = Gcoord(nodemap(localrefineno,:),:);
xl = x(2,1)-x(1,1); % local x length
yh = x(3,2)-x(1,2); % local y height
for ilx = 1:(localx-1)
    xl_coord(ilx+1) = ilx*xl/(localx-1);
end
xl_coord = xl_coord + x(1,1); % local x coordinate
for ily = 1:(localy-1)
    yl_coord(ily+1) = ily*yh/(localy-1);
end
yl_coord = yl_coord + x(1,2); % local y coordinate
lnn = localx*localy; % local number of node
tlnn = 1:1:lnn;
for i = 1:localx
    for j = 1:localy
        n = j+localx*(i-1);
        jseedx = xl_coord(j);
        iseedy = yl_coord(i);
        Lcoord(tlnn(n),:) = [jseedx,iseedy];
    end
end
reference = [1,localx,localx*(localy-1)+1,localx*localy];
Lcoord(reference,:) = [];

Gcoordreduced = Gcoord; %Gcoordreduced(nodemap(localrefineno,:),:)=[];
[ng,mg] = size(Gcoordreduced);[nl,ml] = size(Lcoord);
NewGcoord = zeros((ng+nl),ndim);
NewGcoord(1:ng,:) = Gcoordreduced;
NewGcoord(ng+1:ng+nl,:) = Lcoord;



localnodemap =zeros(ne,eltype);
for inen = 1:(localx-1)
    for layer = 1:(localy-1)
        if layer == 1
            if inen == 1
                %localnodemap(inen+(layer-1)*(localx-1),:) = [inen+(layer-1)*localx, inen+1+(layer-1)*localx, inen+localx+1+(layer-1)*localx, inen+localx+(layer-1)*localx];
                % 1st node of first layer
                localnodemap(1+(layer-1)*(localx-1),:) =  [nodemap(localrefineno,1), nn+1, nn+1+localx-1,nn+1+localx-2]; % 6 17 20 19
            elseif inen == (localx-1)
                % last node of first layer
                localnodemap(localx-1+(layer-1)*(localx-1),:) =  [nn+2,nodemap(localrefineno,2), nn+1+localx+1,nn+1+localx]; % 18 7 22 21
            else
                localnodemap(inen+(layer-1)*(localx-1),:) =  [nn+1, nn+2, nn+1+localx-2+2,nn+1+localx-2+1];
            end
        elseif layer == (localy-1)
            if inen == 1
                localnodemap((localx-1)*(localy-2)+1,:) = [nn+2+(layer-2)*localx+1, nn+2+(layer-2)*localx+2, nn+2+(layer-1)*localx+1, nodemap(localrefineno,4)];
            elseif inen == (localx-1)
                localnodemap((localx-1)*(localy-1),:) = [nn+2+(layer-2)*localx+3, nn+2+(layer-1)*localx, nodemap(localrefineno,3),nn+2+(layer-1)*localx+2];
            else
                localnodemap((localx-1)*(localy-2)+inen,:) = [nn+2+(layer-2)*localx+2, nn+2+(layer-2)*localx+3, nn+2+(layer-1)*localx+2, nn+2+(layer-1)*localx+1];
            end
        else
            localnodemap(inen+(layer-1)*(localx-1),:) = [nn+2+(layer-2)*localx+inen, nn+2+(layer-2)*localx+inen+1, nn+2+(layer-1)*localx+inen+1,nn+2+(layer-1)*localx+inen];
        end
    end
    
end
%nodemap1(n,:) =  [nodemap(localrefineno,i), nn+1, nn+1+localx-1,nn+1+localx-2];
nodemapreduce = nodemap;nodemapreduce(localrefineno,:)=[]; % remove the global node matrix at the local refine area
[ng,mg] = size(nodemapreduce);[nl,ml] = size(localnodemap);
Newnodemap = zeros((ng+nl),eltype);
Newnodemap(1:ng,:) = nodemapreduce;
Newnodemap(ng+1:ng+nl,:) = localnodemap;




