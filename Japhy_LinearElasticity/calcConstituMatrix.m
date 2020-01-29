function D = calcConstituMatrix(E,mu,opt)


if opt == 1                      % plane stress
                
    D = E/(1 - mu*mu) * ...
        [   1      mu        0;
            mu     1         0;
            0      0    (1-mu)/2  ];

elseif opt == 2                  % plane strain
                
    D = E/((1+mu)*(1-2*mu)) * ...
        [  (1-mu)    mu          0;
             mu    (1-mu)        0;
             0       0      (1-2*mu)/2 ];
         
elseif opt == 3                  % axisymmetry
                
    D = E/((1+mu)*(1-2*mu)) * ...
        [  (1-mu)    mu      mu       0;
             mu    (1-mu)    mu       0;
             mu      mu    (1-mu)     0;
             0       0       0    (1-2*mu)/2 ];
end

end   % end of function calcConstituMatrix