% MIRT2D_EXAMPLE2: Non-rigid 2D registration example 2 with Residual Complexity (RC)
% similarity measure, where one of the images is corrupted by nonstationary
% intensity distortion. Try using other similarity measures, e.g. SSD, CC
% or MI to see that they won't work with nonstationary instensity distortions

clear all; close all; clc;
load Sichuan;

im = ImageP(161:440,1:360);
refim_c = ImageMS(41:110,1:90,1:3);
% refim = rgb2gray(refim);
refim_c = imresize(refim_c,size(im));

refim = mean(refim_c,3);

% imshowpair(im,refim,'diff');
tp = refim_c(:,:,1);
refim_c(:,:,1) = refim_c(:,:,3);
refim_c(:,:,3) = tp;

% figure;
% % subplot(1,3,1);
% figure;imshow(refim_c/255);% title('Reference (fixed) image');
% % subplot(1,3,2);
% figure;imshow(im/255); %title('Source image');

 % Main settings
main.similarity='dtv';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI  
main.type = 'TRANSLATION'; %translation
% main.mu=0.01;        % similarity measure parameter (e.g., alpha of RC)
main.TV = TVOP;
main.subdivide=3;       % use 3 hierarchical levels
% main.okno=8;            % mesh window size, the smaller it is the more complex deformations are possible
% main.lambda = 0.005;    % transformation regularization weight, 0 for none
main.single=1;          % show mesh transformation at every iteration
main.alpha = 0.05;
% main.initialX = [1.00 0.10 -5;
%        -0.10 1.00 -0;
%        0 0 1];
% Optimization settings
optim.maxsteps = 200;   % maximum number of iterations at each hierarchical level
optim.fundif = 1e-6;    % tolerance (stopping criterion)
optim.gamma = 1;       % initial optimization step size 
optim.anneal=0.8;       % annealing rate on the optimization step    
 
t0= cputime();
[res1, newim1]=mirt2D_register_rigid(refim,im, main, optim);
t1= cputime()-t0;

main.similarity='rc';

t0= cputime();
[res2, newim2]=mirt2D_register_rigid(refim,im, main, optim);
t2= cputime()-t0;
main.similarity='ssd';

t0= cputime();
[res3, newim3]=mirt2D_register_rigid(refim,im, main, optim);
t3= cputime()-t0;



figure; imshowpair(im, refim,'diff');  title('noregistration');
figure; imshowpair(newim3, refim,'diff'); title('ssd');
figure;imshowpair(newim2, refim,'diff'); title('RC');
figure;imshowpair(newim1, refim,'diff');  title('DTV');

