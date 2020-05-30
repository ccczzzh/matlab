function GQ = calcGaussPoint(nen,NGP)           

if nen == 3               % Triangular element
    if NGP == 1         % One-point quadrature
        GQ.point(1,1) = 1./3;     GQ.point(1,2) = 1./3;
        GQ.weight(1) = 0.5;
    elseif NGP == 3   % Two-point quadrature
        GQ.point(1,1) = 0.5;     GQ.point(1,2) = 0.0;
        GQ.point(2,1) = 0.0;     GQ.point(2,2) = 0.5;
        GQ.point(3,1) = 0.5;     GQ.point(3,2) = 0.5;
        GQ.weight(1) = 1./6;
        GQ.weight(2) = 1./6;
        GQ.weight(3) = 1./6;
    elseif NGP == 4     % Four-point quadrature
        GQ.point(1,1) = 1./3;     GQ.point(1,2) = 1./3;
        GQ.point(2,1) = 0.6;     GQ.point(2,2) = 0.2;
        GQ.point(3,1) = 0.2;     GQ.point(3,2) = 0.6;
        GQ.point(4,1) = 0.2;     GQ.point(4,2) = 0.2;
        GQ.weight(1) = -27./96;
        GQ.weight(2) = 25./96;
        GQ.weight(3) = 25./96;
        GQ.weight(4) = 25./96;
    elseif NGP == 7     % Seven-point quadrature
        GQ.point(1,1) = 1./3;                 GQ.point(1,2) = 1./3;
        GQ.point(2,1) = 0.059715871789770;   GQ.point(2,2) = 0.470142064105115;
        GQ.point(3,1) = 0.470142064105115;   GQ.point(3,2) = 0.059715871789770;
        GQ.point(4,1) = 0.470142064105115;   GQ.point(4,2) = 0.470142064105115;
        GQ.point(5,1) = 0.101286507323456;   GQ.point(5,2) = 0.797426985353087;
        GQ.point(6,1) = 0.101286507323456;   GQ.point(6,2) = 0.101286507323456;
        GQ.point(7,1) = 0.797426985353087;   GQ.point(7,2) = 0.101286507323456;
        GQ.weight(1) = 0.225 / 2;
        GQ.weight(2) = 0.132394152788 / 2;
        GQ.weight(3) = 0.132394152788 / 2;
        GQ.weight(4) = 0.132394152788 / 2;
        GQ.weight(5) = 0.125939180544 / 2;
        GQ.weight(6) = 0.125939180544 / 2;
        GQ.weight(7) = 0.125939180544 / 2;
    else
        disp('ERROR: For triangular elements NGP should be 1, 3, 4 or 7.');
    end
    
elseif nen == 4           % Quadrilateral element
    if NGP == 1         % One-point quadrature
        GQ.point(1,1) = 0.0;  GQ.point(1,2) = 0.0;
        GQ.weight(1) = 4.0;
    elseif NGP == 4     % Four-point quadrature
        GQ.point(1,1) = -sqrt(1/3);   GQ.point(1,2) = -sqrt(1/3);
        GQ.point(2,1) = sqrt(1/3);    GQ.point(2,2) = -sqrt(1/3);
        GQ.point(3,1) = -sqrt(1/3);   GQ.point(3,2) = sqrt(1/3);
        GQ.point(4,1) = sqrt(1/3);    GQ.point(4,2) = sqrt(1/3);
        GQ.weight(1) = 1.0;
        GQ.weight(2) = 1.0;
        GQ.weight(3) = 1.0;
        GQ.weight(4) = 1.0;
    elseif NGP == 9     % Nine-point quadrature
        GQ.point(1,1) = -sqrt(3./5);   GQ.point(1,2) = -sqrt(3./5);
        GQ.point(2,1) = 0.0;          GQ.point(2,2) = -sqrt(3./5);
        GQ.point(3,1) = sqrt(3./5);    GQ.point(3,2) = -sqrt(3./5);
        GQ.point(4,1) = -sqrt(3./5);   GQ.point(4,2) = 0.0;
        GQ.point(5,1) = 0.0;          GQ.point(5,2) = 0.0;
        GQ.point(6,1) = sqrt(3./5);    GQ.point(6,2) = 0.0;
        GQ.point(7,1) = -sqrt(3./5);   GQ.point(7,2) = sqrt(3./5);
        GQ.point(8,1) = 0.0;          GQ.point(8,2) = sqrt(3./5);
        GQ.point(9,1) = sqrt(3./5);    GQ.point(9,2) = sqrt(3./5);
        GQ.weight(1) = 5./9 * 5./9;
        GQ.weight(2) = 8./9 * 5./9;
        GQ.weight(3) = 5./9 * 5./9;
        GQ.weight(4) = 5./9 * 8./9;
        GQ.weight(5) = 8./9 * 8./9;
        GQ.weight(6) = 5./9 * 8./9;
        GQ.weight(7) = 5./9 * 5./9;
        GQ.weight(8) = 8./9 * 5./9;
        GQ.weight(9) = 5./9 * 5./9;
    else
        disp('ERROR: For quadrilateral elements NGP should be 1, 4 or 9.');
    end
end