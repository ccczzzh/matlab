clc;close all;clear all;
prName = input('\nEnter the name of the input file(without the .sol extension): ','s');
inputFile = fopen(strcat(prName,'.sol'), 'r');

fgets(inputFile);  % Read the header line of the file
fgets(inputFile);  % Read the next dummy line.

fgets(inputFile);
for m = 1:18
    dummy = str2num(fgets(inputFile)); 
    x(m,1) = dummy(2);
    uy(m,2) = dummy(5);
    plot(x,uy)
    ylabel('deflection')
    axis([0,4.5,-0.4,0.4])
    hold on;grid on;
    
end
