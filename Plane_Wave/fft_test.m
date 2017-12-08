clear;clc;close all


n = 8;
nx = n;
ny = n;


Vgrid = zeros(nx,ny);
pw = 2;
ix = (ceil((nx-pw)/2)+1):(ceil((nx-pw)/2)+pw)
Vgrid(ceil(nx/2):ceil(nx/2+(pw-1)),(ny/2):(ny/2+(pw-1))) = 1;

Vf = fftshift(fftn(Vgrid)) / (nx*ny);
% Vf = (fftn(Vgrid)) / (nx*ny);

figure(1)
subplot(1,2,1)
b = bar3(Vgrid)
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

subplot(1,2,2)
b = bar3(Vf)
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

% xlim([0,nx])
% ylim([0,ny])
% axis equal