function [lnods_trans] = trans(lnods,nnode,ielem)
%*************************************************************************
% The deriv is defined start from a corner and goes clockwise.
% In this part, just reverse the lnods from clockwise to anti-clockwise to
% make the lnods synic with the derivative, so that the dertiminate will be
% a positive number.
%
% INPUTE ARGUMENTS
% ielem        - The current looping element (according to 'microcell')
% lnods        - Local nodes numbers that consists elements.
% nnode        - Number of nodes concluded in one element.
% 
% OUTPUT 
% lnods_trans  - The transfered lnods.
%
%************************************************************************
lnods_trans = zeros(1,nnode);
if nnode == 3
    lnods_trans(1,1) = lnods(ielem,1);
    lnods_trans(1,2) = lnods(ielem,3);
    lnods_trans(1,3) = lnods(ielem,2);
elseif nnode == 4;
    lnods_trans(1,1) = lnods(ielem,1);
    lnods_trans(1,2) = lnods(ielem,4);
    lnods_trans(1,3) = lnods(ielem,3);
    lnods_trans(1,4) = lnods(ielem,2);
elseif nnode == 6;
    lnods_trans(1,1) = lnods(ielem,1);
    lnods_trans(1,2) = lnods(ielem,6);
    lnods_trans(1,3) = lnods(ielem,5);
    lnods_trans(1,4) = lnods(ielem,4);
    lnods_trans(1,5) = lnods(ielem,3);
    lnods_trans(1,6) = lnods(ielem,2);
elseif nnode == 8;
    lnods_trans(1,1) = lnods(ielem,1);
    lnods_trans(1,2) = lnods(ielem,8);
    lnods_trans(1,3) = lnods(ielem,7);
    lnods_trans(1,4) = lnods(ielem,6);
    lnods_trans(1,5) = lnods(ielem,5);
    lnods_trans(1,6) = lnods(ielem,4);
    lnods_trans(1,7) = lnods(ielem,3);
    lnods_trans(1,8) = lnods(ielem,2);
elseif nnode == 9;
    lnods_trans(1,1) = lnods(ielem,1);
    lnods_trans(1,2) = lnods(ielem,8);
    lnods_trans(1,3) = lnods(ielem,7);
    lnods_trans(1,4) = lnods(ielem,6);
    lnods_trans(1,5) = lnods(ielem,5);
    lnods_trans(1,6) = lnods(ielem,4);
    lnods_trans(1,7) = lnods(ielem,3);
    lnods_trans(1,8) = lnods(ielem,2);
    lnods_trans(1,9) = lnods(ielem,9);
end
end