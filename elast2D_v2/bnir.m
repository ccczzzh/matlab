% first of, need to find the positioin of all the integration points
n = 6; % here need to type in the number of integration points
pos = zeros(n,1);
for i = 1:n
    pos(i) = -1 + (1/n) + (2/n)*(i-1);
end
syms x 
A = zeros(n,1);
for i = 1:n
    A(i) = int(x^(i-1),-1,1);
end
M = zeros(n,n);
for i = 1:n
    for j = 1:n
        M(i,j) = pos(j)^(i-1);
    end
end
W = M\A;
M = W*W';
N = n^2;
posip = zeros(2,N);
for i = 1:n
    for j = 1:n
        posip(1,(i-1)*n+j) = pos(i);
        posip(2,(i-1)*n+j) = pos(j);
    end
end
weiip = zeros(N,1);
for i = 1:n
    for j = 1:n
        weiip((i-1)*n+j) = M(i,j);
    end
end

