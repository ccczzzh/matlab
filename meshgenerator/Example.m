clc;clear all;close all;

Lx = input('Input length of x: ','s');Lx = str2num(Lx);
Ly = input('Input length of y: ','s');Ly = str2num(Ly);

Nx = input('Input number of element on x: ','s');Nx = str2num(Nx);
Ny = input('Input number of element on y: ','s');Ny = str2num(Ny);

FEMN = input('Enter the number of nodes for Mesh: ','s'); FEMN = str2num(FEMN);
tic

    
if FEMN == 3
    NE = 2*Ny*Nx;
    [coords nT nNodes]=Mesh3nodes(Lx,Ly,Nx,Ny);
    
    disp(['Number of Nodes =  ', num2str(nNodes)])
    disp('Connectivity Table')
    disp(nT)
    
    z = 1;
    for i = 1:NE
        figure(1); patch('Vertices',coords(z:z+2,:),'Faces',[1,2,3],'FaceColor'...
            ,'none','EdgeColor','g')
        hold on;
        z = z+3;
    end
    
    figure(1),scatter(coords(:,1),coords(:,2),'MarkerFaceColor','r')
    hold off;
    
elseif FEMN == 4
    NE = Ny*Nx;
    [coords nT nNodes]=Mesh4nodes(Lx,Ly,Nx,Ny);
    
    disp(['Number of Nodes =  ', num2str(nNodes)])
    disp('Connectivity Table')
    disp(nT)
    
    z = 1;
    for i = 1:NE
        figure(1); patch('Vertices',coords(z:z+3,:),'Faces',[1,2,4,3],'FaceColor'...
            ,'none','EdgeColor','g')
        hold on;
        z = z+4;
    end
    
    figure(1),scatter(coords(:,1),coords(:,2),'MarkerFaceColor','r')
    hold off;

end
toc