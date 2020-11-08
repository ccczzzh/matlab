function [BW]=bandwidth(ne,nen,node,ndim)
% Calculate the Bandwidth
BW = 0;
for i = 1:ne
   % for j = 1:nen
        MABW = node(i,1);
        MIBW = node(i,1);
        for k = 2:nen
            temp = node(i,k);
            if MABW < temp
                MABW = temp;
            end
            if MIBW > temp
                MIBW = temp;
            end
            iBW = (MABW-MIBW+1)*ndim;
            if BW<iBW
                BW = iBW
            end
        end
 %   end
end
