clc;clear all;close all;
% ========================================
% simualtion free fall ball
% ========================================
% free fall eqution yt = v0*t+0.5*(-g)*t^2
% free fakk velocity vt = v0 +g*t
% unit: m, s,
% time difference
n = 0.02;
% bounce ground back ratio
kn = -0.8;
% initial vertical location
yv0 = 10;
%initial horizontal location
yh0 = 0;
% initial vertical velocity
vv0 = 35;
% initial horizontal velocity
vh0 = 4;
% gravity
g = -9.81;
yvt = zeros(1000,1);
yht = zeros(1000,1);
tm = zeros(1000,1);
vt = zeros(1000,1);
for t = 0:1:999
%     yht(1) = yh0;
%     yht(0) = 0;
    yvt(t+1) = yv0+vv0*(n)+0.5*g*((n)^2);
    yv0 = yvt(t+1);
    yht(t+1) = vh0*n*(t+1);
    yh0 = yht(t+1);
    tm(t+1) = (t);
    vt(t+1) = vv0 + g*n;
    vv0 = vt(t+1);
       
    if yv0 <= 0
        vv0 = kn*vt(t+1);
        %break
    end
end
figure(1); title('Bounce Ball Simulator');
tiledlayout(2,1)
ax1 = nexttile;
plot(ax1,yht,yvt,'x');
xlabel('time');ylabel('displacement');
ax2 = nexttile;
plot(ax2,tm,vt,'o');
xlabel('time');ylabel('velocity');

