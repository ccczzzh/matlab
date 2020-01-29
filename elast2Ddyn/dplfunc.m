%***********************************************************************
% function dplfun
%
% DESCRIPTION
% Derivative of the piecewise linear function defined by a set of npoint
% pairs {X,F(X)} stored in the matrix xfx. (dim. 2*npont).
%
% HISTORY
% E A de Souza Neto/M Partovi,               April 2004: Initial coding
%***********************************************************************
%
function [dplfun] = dplfunc(x,npoint,xfx)
%
%***********************************************************************
if x >= xfx(1,npoint)
% x >= x(npoint) --> f(x) = f(x(npoint)) -->df/dx = 0
   dplfun = 0;
else
    for i=1:npoint
        if x >= xfx(1,i)
           continue 
        elseif i == 1
%       ---x< x1 --> f(x)=f(x1)--> df(x)/dx = 0
            dplfun = 0;
            break
        else 
%       ---x(i-1) <= x < x(i) 
         dplfun = (xfx(2,i)-xfx(2,i-1))/(xfx(1,i)-xfx(1,i-1));
        end
    end
end
