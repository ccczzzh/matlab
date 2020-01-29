function [K,Fn] = assembly(K,Fn,ke,lm0,lm,u,eldof,dbcVal)


for ii = 1:eldof
    P = lm(ii);    % P is the global position of ii
    if P > 0
        for jj = 1:eldof
            Q = lm(jj);  % Q is the global postion of jj
            if Q > 0
                K(P,Q) = K(P,Q) + ke(ii,jj);
            else
                Q0 = lm0(jj);
                Fn(P) = Fn(P) - ke(ii,jj) * (dbcVal(dbcVal==Q0,2) - u(Q0));
            end
        end
    end
end


end  % end of function assemble