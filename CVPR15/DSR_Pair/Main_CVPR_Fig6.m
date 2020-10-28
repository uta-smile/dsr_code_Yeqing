clear all; close all; clc;
load mirt2D_data2.mat;  
load T;


im=mirt2D_transform(refim, res1);
im = add_distor(im,6);

 % Main settings
main.similarity='dtv';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI  
% main.mu=0.01;        % similarity measure parameter (e.g., alpha of RC)
main.TV = TVOP;
main.subdivide=3;       % use 3 hierarchical levels
main.okno=8;            % mesh window size, the smaller it is the more complex deformations are possible
main.lambda = 0.005;    % transformation regularization weight, 0 for none
main.single=0;          % show mesh transformation at every iteration

% Optimization settings
optim.maxsteps = 200;   % maximum number of iterations at each hierarchical level
optim.fundif = 1e-5;    % tolerance (stopping criterion)
optim.gamma = 1;       % initial optimization step size 
optim.anneal=0.8;       % annealing rate on the optimization step    
 
tic;
[res1, newim1]=mirt2D_register(refim,im, main, optim);
tdtv = toc;


% Main settings
main.similarity='rc';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI  
main.alpha=0.05;        % similarity measure parameter (e.g., alpha of RC)

tic;
[res2, newim2]=mirt2D_register(refim,im, main, optim);
trc = toc;

main.similarity='ssd';
tic;
[res3, newim3]=mirt2D_register(refim,im, main, optim);
tssd = toc;

fprintf('SSD Running time %f\n', tssd);
fprintf('RC Running time %f\n', trc);
fprintf('DTV Running time %f\n', tdtv); 

figure;
subplot(2,4,1);imshow(refim); title('Reference (fixed) image');
subplot(2,4,5);imshow(im);    title('Source (float) image');
subplot(2,4,4);imshow(newim1); title('Registered by dTV');
subplot(2,4,8);mirt2D_meshplot(res1.X(:,:,1),res1.X(:,:,2));
subplot(2,4,3);imshow(newim2); title('Registered by RC');
subplot(2,4,7);mirt2D_meshplot(res2.X(:,:,1),res2.X(:,:,2));
subplot(2,4,2);imshow(newim3); title('Registered by SSD');
subplot(2,4,6);mirt2D_meshplot(res3.X(:,:,1),res3.X(:,:,2));

% figure;mirt2D_meshplot(main.X(:,:,1),main.X(:,:,2));

gradref = gradientimage(refim);
gradim = gradientimage(im);
grad1 = gradientimage(newim1);
grad2 = gradientimage(newim2);
figure;
subplot(2,4,5);imshow(gradim,[]);    title('Source (float) image');
subplot(2,4,1);imshow(gradref,[]); title('Reference (fixed) image');
subplot(2,4,4);imshow(grad1,[]); title('Registered by dTV');
subplot(2,4,3);imshow(grad2,[]); title('Registered by RC');