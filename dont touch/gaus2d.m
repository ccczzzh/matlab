function [posgp,weigp] = gaus2d(domain,ngaus)
% Sets Gaussian quadrature constants for 2-D domains
%***********************************************************************
% DESCRIPTION
% Given the type of domain and the required number of integration
% points, this routine sets the sampling point positions and the
% corresponding weights following standard Gauss quadrature rules
% for 2-D domains (quadrilateral and triangular domains).
%
% INPUT
% domain  - domain type flag.(string)
%              'QUA' for quadrilateral domains
%              'TRI' for triangular domains
% ngaus    - number of integration points.
%
% OUTPUT
% posgp    - array containing the positions of the integration points.
% weigp    - array containing the weights at the integration points.
%
% HISTORY
% M Partovi, October 2004: Initial coding (following similar routine
%                          with same name in program HYPLAS)
%***********************************************************************

%***********************************************************************
% SET SAMPLING POINTS POSITIONS AND WEIGHTS FOR GAUSSIAN NUMERICAL
% INTEGRATION RULES IN 2-D
%***********************************************************************
%
if strcmp(domain,'QUA')   
%
% Integration over quadrilateral domain with vertices
%
% (-1,1)   _ _ _ _ _ _ _ _ _ _ _  (1,1)
%         |          |          |
%         |          |          |
%         |          |          |
%         |_ _ _ _ _ |_ _ _ _ _ |
%         |          |          |  
%         |          |          |
%         |          |          |
% (-1,-1) |_ _ _ _ _ |_ _ _ _ _ | (1,-1)
%       
%***********************************************************************
   if ngaus == 1
      posgp(1,1) = 0;
      posgp(2,1) = 0;
      weigp(1)   = 4;
   elseif ngaus == 4
          posgp(1,1) = -0.577350269189626;
          posgp(2,1) = -0.577350269189626;
          weigp(1)   = 1.0;
          posgp(1,2) = -0.577350269189626;
          posgp(2,2) = +0.577350269189626;
          weigp(2)   = 1.0;
          posgp(1,3) = +0.577350269189626;
          posgp(2,3) = -0.577350269189626;
          weigp(3)   = 1.0;
          posgp(1,4) = +0.577350269189626;
          posgp(2,4) = +0.577350269189626;
          weigp(4)   = 1.0;
   elseif ngaus == 5
          posgp(1,1) = -0.774596669241483;
          posgp(2,1) = -0.774596669241483;
          weigp(1)   = 0.5555555555555556;
          posgp(1,2) = -0.774596669241483;
          posgp(2,2) = +0.774596669241483;
          weigp(2)   = 0.5555555555555556;
          posgp(1,3) = +0.774596669241483;
          posgp(2,3) = -0.774596669241483;
          weigp(3)   = 0.5555555555555556;
          posgp(1,4) = +0.774596669241483;
          posgp(2,4) = +0.774596669241483;
          weigp(4)   = 0.5555555555555556;
          posgp(1,5) = +0.0;
          posgp(2,5) = +0.0;
          weigp(5)   = 1.7777777777777778;
   elseif ngaus == 9
          posgp(1,1) = -0.774596669241483;
          posgp(2,1) = -0.774596669241483;
          weigp(1)   = +0.308641975308643;
          posgp(1,2) = -0.774596669241483;
          posgp(2,2) = +0.0;
          weigp(2)   = +0.493827160493828;
          posgp(1,3) = -0.774596669241483;
          posgp(2,3) = +0.774596669241483;
          weigp(3)   = +0.308641975308643;
          posgp(1,4) = +0.0;
          posgp(2,4) = -0.774596669241483;
          weigp(4)   = +0.493827160493828;
          posgp(1,5) = +0.0;
          posgp(2,5) = +0.0;
          weigp(5)   = +0.790123456790124;
          posgp(1,6) = +0.0;
          posgp(2,6) = +0.774596669241483;
          weigp(6)   = +0.493827160493828;
          posgp(1,7) = +0.774596669241483;
          posgp(2,7) = -0.774596669241483;
          weigp(7)   = +0.308641975308643;
          posgp(1,8) = +0.774596669241483;
          posgp(2,8) = +0.0;
          weigp(8)   = +0.493827160493828;
          posgp(1,9) = +0.774596669241483;
          posgp(2,9) = +0.774596669241483;
          weigp(9)   = +0.308641975308643;
   end
%          
elseif strcmp(domain,'TRI') 
%
% Integration over triangular domain with vertices 
%
%       (0,1)
%         |\
%         | \
%         |  \
%         |   \
%         |    \
%         |     \
%         |_ _ _ \
%   (0,0)          (1,0)
%
%***********************************************************************
   if ngaus == 1
      posgp(1,1) = 0.333333333333333;
      posgp(2,1) = 0.333333333333333;
      weigp(1)   = 0.5;
   elseif ngaus == 3
          posgp(1,1) = 0.166666666666667;
          posgp(2,1) = 0.166666666666667;
          weigp(1)   = 0.166666666666667;
          posgp(1,2) = 0.666666666666667;
          posgp(2,2) = 0.166666666666667;
          weigp(2)   = 0.166666666666667;
          posgp(1,3) = 0.166666666666667;
          posgp(2,3) = 0.666666666666667;
          weigp(3)   = 0.166666666666667;
   end
end