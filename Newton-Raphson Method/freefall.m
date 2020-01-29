
tic;
step=1000;
%dt=5/step;
dt=15/step;

y_out=[];
t_out=[];

%initialize time loop
time=0;
yn=0;

for i=1:step
time=time+dt

%integrate ( solve for y (n+1) )
y = BENR(yn,dt);        %3 newton raphson

%store results
t_out(i)=time;
y_out(i)=y;
%y_out(i)=(y-yn)/dt;

% update time step
yn=y;

end;

hold all
plot(t_out,y_out);
toc;

    





