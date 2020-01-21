clc;close all;clear all;

A = [1 2 3;4 5 6;7 8 9]

B = zeros(9,1)

for t = 0:1:8
    B(t+1) = 4*t-6+t/5
end
C = rand(3)

