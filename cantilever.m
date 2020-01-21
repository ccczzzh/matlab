clc;clear all;close all;

% cantilever property
% length
L = 10; dl = L/10;
% Young's Modulus
E = 200e6;
% Section shape: 1 = square; 2 = circular; 3 = hollow circular.
section = 1;
a = 10;
I = areaMomentofInertia(10,0 ,0 , 1);

% load 
F = 10;
%working length
y = 0:dl:L;
n = length(y);
dB = zeros(11,1); dL = zeros(11,1);dM = zeros(11,1);
for l = 0:dl:L
   dB(l+1) = -(F*l^3)/(3*E*I);
   dL (l+1) = l;
   dM(l+1) = -F*l;
end
figure(1);tiledlayout(2,1);ax1 = nexttile;
plot(ax1,dL,dB);
xlabel('Location');ylabel('Bending displacement');
ax2 = nexttile;
plot(ax2,dM,dL);
xlabel('Location');ylabel('Bending Moment');
