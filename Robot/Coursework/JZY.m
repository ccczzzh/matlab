%% e 2 Link
figure;
% Find Link DH parameters, the coordinate set up can be found on report
L(1)=Link([0,1,0,pi/2]);
L(2)=Link([0,0,1,-pi/2]);
L(1).offset=0;
L(2).offset=0;
RR=SerialLink(L);
%name the robot
RR.name = 'LOL';
% Initial pose [pi/2, pi/4].
RR.plot([pi/2,pi/4]);
%% f Draw workspace
% assume the working range of the robot
figure;
L(1).qlim=[-pi,pi];
L(2).qlim=[-pi/4,pi/4];
% label axes and unit
xlabel('X (meters)'); ylabel('Y (meters)'); zlabel('Z(meters)');
title('Workspace of RR'); hold on; N=40;
% Plot the robot with initial pose [0 0]
for i= 1:N+1
for j= 1:N+1
% These value should be setting basing on the link limit so that
% the link can rotate from the minimum to the maximum value however
% this question has not stated working range therefore follow the assumption
TR= RR.fkine([-pi+2*pi*(i-1)/N -pi/4+pi/2*(j-1)/N]);
% Using the coordinate of the end effector (column vector) to creat a 3D matrix
SURF(i,j,:) = TR(:,4); % TR(:,4), translation vector
end;
end;
% Using the 3D matrix to plot the location of the end effector has reached
%to form a surface of the reachable working space
surf(SURF(:,:,1), SURF(:,:,2), SURF(:,:,3));
RR.plot([pi/2,pi/4]);