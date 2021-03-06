%***********************************************************************
% function plfun
%
% DESCRIPTION
% Piecewise linear function defined by a set of npoint
% pairs {X,F(X)} stored in the matrix xfx. (dim. 2*npont).
%
% HISTORY
% E A de Souza Neto/M Partovi,  April 2004: Initial coding
% EAdSN                         March 2015: Bug fix
%***********************************************************************
%
function [plfun] = plfunc(x,npoint,xfx)
%
%***********************************************************************
if x >= xfx(1,npoint)
% ---  x >= x(npoint) --> f(x) = f(x(npoint))  ---
   plfun = xfx(2,npoint);
else
   for i=1:npoint    
       if x >= xfx(1,i)
          continue        
       elseif i == 1
%         -- x < x1 --> f(x)=f(x1) ---
          plfun = xfx(2,1);
          break
       else
%         -- x(i-1) <= x < x(i) ---
          plfun = xfx(2,i-1) + (x-xfx(1,i-1))*...
          (xfx(2,i)-xfx(2,i-1))/(xfx(1,i)-xfx(1,i-1));
          break
       end
   end
end