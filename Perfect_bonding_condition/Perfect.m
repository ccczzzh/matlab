%PerfectBondingCondition
D1 = xlsread('perfect.xlsx','A:A');
F1 = xlsread('perfect.xlsx','B:B');
hold on
plot(D1,F1)

legend('Perfet Bonding','Location','northeast')