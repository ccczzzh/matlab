% clc;clear all;close all;
%read the Excel File
CE = xlsread('Concrete_E.xlsx');
D1 = CE(:,1);F1 = CE(:,2);
D2 = CE(:,3);F2 = CE(:,4);
D3 = CE(:,5);F3 = CE(:,6);
D4 = CE(:,7);F4 = CE(:,8);
figure;
plot(D1,F1);
hold on;
plot(D2,F2);plot(D3,F3);plot(D4,F4);
title('Concrete Behavior')
xlabel('Displacement')
ylabel('Load')
% xlim([0 0.017])
% ylim([-0.2 0.01])
legend('CE D=1e10','CE D=1','cE D=1e3','CE D=1e7','Location','northeast')

CE = xlsread('Steel_E.xlsx');
SD1 = CE(:,1);SF1 = CE(:,2);
SD2 = CE(:,3);SF2 = CE(:,4);
SD3 = CE(:,5);SF3 = CE(:,6);
SD4 = CE(:,7);SF4 = CE(:,8);
figure;
plot(SD1,SF1);
hold on;
plot(SD2,SF2);plot(SD3,SF3);plot(SD4,SF4);
title('Steel Behavior')
xlabel('Displacement')
ylabel('Load')
% xlim([0 0.017])
% ylim([-0.2 0.01])
legend('All D=1e10','All D=1','All D=1e3','All D=1e7','Location','northeast')
