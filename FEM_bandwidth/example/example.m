x = [0 1 1 0 0];
y = [0 0 1 1 0];
stress = [1 2 3 4 1];
fill(x,y,stress)
shading interp;
colorbar;
axis equal
%%
t = linspace(0,pi)';
x = cos(t);
y = sin(t);
color = jet(7);
figure
hold on
for kk = 1:7
  plot((2+kk)*x,(2+kk)*y,'Color',color(kk,:),'LineWidth',20)
end
ylim([0 15])
%%
[xGrid,yGrid] = meshgrid(-10:0.1:10,0:0.1:10);
zGrid = sqrt(xGrid.^2 + yGrid.^2);
figure
surf(xGrid,yGrid,zGrid,...
  'EdgeColor',  'none',...
  'FaceAlpha',  0.5)
colormap(jet)
set(gca,'CLim',[4 11])
zlim([5 10])
view(2)
%%
i = [1 2 3 4 5 6 7 8 9 10];
for j= 1:10
    if i(j) >7
        i(j)
        return
    end
end
