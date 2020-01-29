function [deriv] = SFandD(s,t,nnode)
% Set the shape function and the corresponding derivatives based on the
% given number of nodes
% the s(xi) and t(eta) will read in the function stiffness
%
% INPUT ARGUMENTS
% t,s     - corresponding number following the gauss quadrature rule
% nnode   - number of nodes concluded in one element
%
% OUTPUT
% deriv   - two by nnode array of derivatives
%
%***********************************************************************
    ss = s*s;
    tt = t*t;
    st = s*t;
    sst = s*s*t;
    stt = s*t*t;
    sstt = s*s*t*t;
if nnode == 4;
    % isoparametric shape function values at gauss point
    % nodes are logged clockwise start from the bottom left corner
    shape(1) = (1-t-s+st)*0.25;
    shape(2) = (1-t+s-st)*0.25;
    shape(3) = (1+t+s+st)*0.25;
    shape(4) = (1+t-s-st)*0.25;
    % and isoparametric derivatives
    deriv(1,1) = (-1+t)*0.25;
    deriv(1,2) = (+1-t)*0.25;
    deriv(1,3) = (+1+t)*0.25;
    deriv(1,4) = (-1-t)*0.25;
    deriv(2,1) = (-1+s)*0.25;
    deriv(2,2) = (-1-s)*0.25;
    deriv(2,3) = (+1+s)*0.25;
    deriv(2,4) = (+1-s)*0.25;
elseif nnode == 6;
    % nodes are logged clockwise start any corner
    shape(1) = 2*ss-s;
    shape(2) = 4*st;
    shape(3) = 2*tt-t;
    shape(4) = 4*t-4*tt-4*st;
    shape(5) = 1-3*s-3*t+2*ss+2*tt+4*st;
    shape(6) = 4*s-4*st-4*ss;
    % and isoparametric derivatives
    deriv(1,1) = 4*s-1;
    deriv(2,1) = 0;
    deriv(1,2) = 4*t;
    deriv(2,2) = 4*s;
    deriv(1,3) = 0;
    deriv(2,3) = 4*t-1;
    deriv(1,4) = -4*t;
    deriv(2,4) = 4-8*t-4*s;
    deriv(1,5) = -3+4*s+4*t;
    deriv(2,5) = -3+4*t+4*s;
    deriv(1,6) = 4-4*t-8*s;
    deriv(2,6) = -4*s;
elseif nnode == 8;
    % nodes are logged clockwise start from the bottom left corner
    shape(1) = (-1+st+ss+tt-sst-stt)*0.25;
    shape(2) = (1-t-ss+sst)*0.5;
    shape(3) = (-1-st+ss+tt-sst+stt)*0.25;
    shape(4) = (1+s-tt-stt)*0.5;
    shape(5) = (-1+st+ss+tt+sst+stt)*0.25;
    shape(6) = (1+t-ss-sst)*0.5;
    shape(7) = (-1-st+ss+tt+sst-stt)*0.25;
    shape(8) = (1-s-tt+stt)*0.5;
    % and isoparametric derivatives
    deriv(1,1) = (t+s*2-st*2-tt)*0.25;
    deriv(1,2) = -s+st;
    deriv(1,3) = (-t+s*2-st*2+tt)*0.25;
    deriv(1,4) = (1-tt)*0.5;
    deriv(1,5) = (t+s*2+st*2+tt)*0.25;
    deriv(1,6) = -s-st;
    deriv(1,7) = (-t+s*2+st*2-tt)*0.25;
    deriv(1,8) = (-1+tt)*0.5;
    deriv(2,1) = (s+t*2-ss-st*2)*0.25;
    deriv(2,2) = (-1+ss)*0.5;
    deriv(2,3) = (-s+t*2-ss+st*2)*0.25;
    deriv(2,4) = -t-st;
    deriv(2,5) = (s+t*2+ss+st*2)*0.25;
    deriv(2,6) = (1-ss)*0.5;
    deriv(2,7) = (-s+t*2+ss-st*2)*0.25;
    deriv(2,8) = -t+st;
elseif nnode == 9;
    % nodes are logged clockwise start from the bottom left corner
    shape(1) = (st-stt-sst+sstt)*0.25;
    shape(2) = (-t+tt+sst-sstt)*0.5;
    shape(3) = (-st+stt-sst+sstt)*0.25;
    shape(4) = (s+ss-stt-sstt)*0.5;
    shape(5) = (st+sst+stt+sstt)*0.25;
    shape(6) = (t-sst+tt-sstt)*0.5;
    shape(7) = (-st-stt+sst+sstt)*0.25;
    shape(8) = (-s+stt+ss-sstt)*0.5;
    shape(9) = 1-ss-tt+sstt;
    % shape functions derivatives
    deriv(1,1) = (t-tt-st*2+stt*2)*0.25;
    deriv(1,2) = st-stt;
    deriv(1,3) = (-t+tt-st*2+stt*2)*0.25;
    deriv(1,4) = (1+s*2-tt-stt*2)*0.5;
    deriv(1,5) = (t+st*2+tt+stt*2)*0.25;
    deriv(1,6) = -st-stt;
    deriv(1,7) = (-t-tt+st*2+stt*2)*0.25;
    deriv(1,8) = (-1+tt+s*2-stt*2)*0.5;
    deriv(1,9) = -s*2+stt*2;
    deriv(2,1) = (s-st*2-ss+sst*2)*0.25;
    deriv(2,2) = (-1+t*2+ss-sst*2)*0.5;
    deriv(2,3) = (-s+st*2-ss+sst*2)*0.25;
    deriv(2,4) = -st-sst;
    deriv(2,5) = (s+ss+st*2+sst*2)*0.25;
    deriv(2,6) = (1-ss+t*2-sst*2)*0.5;
    deriv(2,7) = (-s-st*2+ss+sst*2)*0.25;
    deriv(2,8) = st-sst;
    deriv(2,9) = -t*2+sst*2;
end