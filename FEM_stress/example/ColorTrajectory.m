clc;close all;clear all;
t = linspace(0,2*pi,100)
x = cos(t)
y = sin(t)
z = t
c = cos(t).^2

colormap(hsv)
patch(x,y,z,c,'FaceColor','none','EdgeColor','interp')
colorbar
view(3)