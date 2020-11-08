function [UXnew,UYnew] = solveandresult(nq,dn,K,F,U,nn)
move = setdiff(1:nq,dn);
Kreduced = K(move,move);
Freduced = F(move);
Ureduced = Kreduced\Freduced;
U(move) = Ureduced;
F = K*U;
Unew = reshape(U,[2,nn]);
Unew = Unew.';
% Deformed X Y Coordinate
UXnew = Unew(:,1);UYnew = Unew(:,2);