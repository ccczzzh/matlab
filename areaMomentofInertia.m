% calculate are moment of inertia
function I=amoi(a,R,r,section)
if section == 1
    % cross section is square
    I = (a^4)/12;
elseif section == 2
    % cross section is circular 
    I = (pi()*r^4)/4;
elseif section == 3
    % cross section is hollow circular
    I = (pi()*(R^4-r^4))/64;
end

end