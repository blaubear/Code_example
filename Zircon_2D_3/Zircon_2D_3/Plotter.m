clc
nx = 512;
ny = nx;
figure;
fileID = fopen('2.txt','r');
formatSpec = '%f';
sizeA = [nx nx];
A = fscanf(fileID,formatSpec,sizeA);

Lx = 1;
Ly = Lx;

c = colorbar, c.Limits=[0 1];
dx = Lx/(nx-1);
dy = Ly/(ny-1);
[X,Y] = ndgrid(0:dx:Lx,0:dy:Ly);
pcolor(X,Y,A),shading flat,axis image;c = colorbar, c.Limits=[0 500];drawnow
xlabel('x');
ylabel('y');