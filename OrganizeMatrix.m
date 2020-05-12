clc;close all;clear all;
tic;
inp = rand(10,1);
[a,b] = size(inp);
out = inp;
for i =1:a
    for j = i+1:a
        if out(i,1)>out(j,1)
            out([i j]) = out([j i]);
          
 
        end
    end

end
toc;
a = 1:1:a;
figure;
bar(a,inp);
figure;
bar(a,out);