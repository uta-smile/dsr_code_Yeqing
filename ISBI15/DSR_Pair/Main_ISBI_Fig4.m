
clear all; close all; clc;

fixed = rgb2gray(imread('iris1.jpg'));
moving = rgb2gray(imread('iris2.jpg'));

refim_dis = double(fixed);
refim_dis = imresize(refim_dis,[200,250]);
refim_dis = refim_dis/max(refim_dis(:));

im = double(moving);
im = imresize(im,[200,250]);
im = im/max(im(:));

% figure;
% subplot(2,3,1);imshow(refim_dis); title('Reference (fixed) image');
% subplot(2,3,2);imshow(im); title('Source image');
% 
% figure, imshowpair(refim_dis,im, 'montage');


 % Main settings
main.similarity='dtv';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI  
% main.mu=0.01;        % similarity measure parameter (e.g., alpha of RC)
main.TV = TVOP;
main.subdivide=4;       % use 3 hierarchical levels
main.okno=8;            % mesh window size, the smaller it is the more complex deformations are possible
main.lambda = 0.005;    % transformation regularization weight, 0 for none
main.single=0;          % show mesh transformation at every iteration

% Optimization settings
optim.maxsteps = 200;   % maximum number of iterations at each hierarchical level
optim.fundif = 1e-7;    % tolerance (stopping criterion)
optim.gamma = 1;       % initial optimization step size 
optim.anneal=0.8;       % annealing rate on the optimization step    
 
t0=cputime();
[res1, newim1]=mirt2D_register(refim_dis,im, main, optim);
t1=cputime()-t0;

% Main settings
main.similarity='rc';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI  
main.alpha=0.05;        % similarity measure parameter (e.g., alpha of RC)

 t0=cputime();
[res2, newim2]=mirt2D_register(refim_dis,im, main, optim);
t2=cputime()-t0;

figure, imshowpair(im, refim_dis)
title('Unregistered')

figure, imshowpair(newim2, refim_dis);
title('RC')

figure, imshowpair(newim1, refim_dis);
title('DTV')


