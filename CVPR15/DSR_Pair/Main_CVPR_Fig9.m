% MIRT2D_EXAMPLE2: Non-rigid 2D registration example 2 with Residual Complexity (RC)
% similarity measure, where one of the images is corrupted by nonstationary
% intensity distortion. Try using other similarity measures, e.g. SSD, CC
% or MI to see that they won't work with nonstationary instensity distortions

clear all; close all; clc;

im = double(imread('westconcordaerial.png'));
refim = double(imread('westconcordorthophoto.png'));

im = im/255;
refim = refim/255;
im_c = im;
im = rgb2gray(im);
% centerFixed = size(fixedVolume)/2;
% centerMoving = size(movingVolume)/2;
% figure, title('Unregistered Axial slice');
% imshowpair(movingVolume(:,:,centerMoving(3)), );

% im = movingVolume(:,:,centerMoving(3));
% refim = fixedVolume(:,:,centerFixed(3));
refimlarge = zeros(size(im));
refimlarge(1:366,1:364) = refim;
refim = refimlarge;

% imshowpair(im,refim);

% 
% refim = rgb2gray(refim_c);
% [m,n] = size(refim);
% % refim = double(refim)/255;
% 
% im_c = imresize(im_c,[m,n]);
% im = rgb2gray(im_c);

% [m,n] = size(im);
% im = double(im)/255;

% refim = imresize(refim,[128,128]);
% I = zeros(150,150);
% I(12:end-11,12:end-11) = refim;
% refim = I;
% Tim = [1.02 0.03 -2
%        -0.03 1.03 -3];

% im = warpAffine2(refim,Tim);

% im = im_int(:,:,1);
% refim_dis = add_distor(refim,2,0);
% refim_dis = add_distor(refim_dis,2,1);
% 
% figure;
% subplot(2,1,1);imshow(refim); title('Reference (fixed) image');
% subplot(2,1,2);imshow(im); title('Source image');


 % Main settings
main.similarity='dtv';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI  
main.type = 'AFFINE'; %translation
% main.mu=0.01;        % similarity measure parameter (e.g., alpha of RC)
main.TV = TVOP;
main.subdivide=3;       % use 3 hierarchical levels
% main.okno=8;            % mesh window size, the smaller it is the more complex deformations are possible
% main.lambda = 0.005;    % transformation regularization weight, 0 for none
main.single=0;          % show mesh transformation at every iteration
main.alpha = 0.05;
main.initialX = [1.00 0.10 -5;
       -0.10 1.00 -0;
       0 0 1];
% Optimization settings
optim.maxsteps = 200;   % maximum number of iterations at each hierarchical level
optim.fundif = 1e-6;    % tolerance (stopping criterion)
optim.gamma = 1;       % initial optimization step size 
optim.anneal=0.8;       % annealing rate on the optimization step    
 
opts.type = 'AFFINE'; %translation
% opts.TV = TVOP;
opts.levels=3;       % use 3 hierarchical levels
opts.display=1;          % show mesh transformation at every iteration
opts.maxsteps = 200;   % maximum number of iterations at each hierarchical level
opts.fundif = 1e-6;    % tolerance (stopping criterion)
opts.gamma = 1;       % initial optimization step size 
opts.anneal=0.8;       % annealing rate on the optimization step    
opts.dimen = 2; % 2D
 opts.initialX = [1.00 0.10 -5;
       -0.10 1.00 -0;
       0 0 1];



t0= cputime();
[res1, newim1]=Register_DTV(refim,im, opts);
t1= cputime()-t0;

main.similarity='rc';

t0= cputime();
[res2, newim2]=mirt2D_register_rigid(refim,im, main, optim);
t2= cputime()-t0;
main.similarity='ssd';

t0= cputime();
[res3, newim3]=mirt2D_register_rigid(refim,im, main, optim);
t3= cputime()-t0;

%%

% figure,subplot(2,3,1); imshow( refim_c);  title('ref');
% subplot(2,3,2); imshow(im_c);  title('source');
% subplot(2,3,3); imshow(warpAffine2(im_c,res1.X)); title('DTV');
% subplot(2,3,4); imshow(warpAffine2(im_c,res2.X)); title('RC');
% subplot(2,3,5); imshow(warpAffine2(im_c,res3.X)); title('SSD');
% % % Main settings
% % main.similarity='rc';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI  
% % main.alpha=0.05;        % similarity measure parameter (e.g., alpha of RC)
% figure,subplot(2,2,1); imshowpair(im, refim);  title('UnRegistered');
% subplot(2,2,4); imshowpair(newim3, refim); title('SSD');
% subplot(2,2,3); imshowpair(newim2, refim); title('RC');
% subplot(2,2,2); imshowpair(newim1, refim); title('DTV');


% figure; imshow(refim_c);  
% figure;imshow(im_c); 
% figure; imshow(warpAffine2(im,res3.X));
% figure; imshow(warpAffine2(im,res2.X)); 
% figure; imshow(warpAffine2(im,res1.X)); 

%--------------------------------
% im = flip(im, 1);
% refim = flip(refim, 1);
% newim1 = flip(newim1, 1);
% newim2 = flip(newim2, 1);
% newim3 = flip(newim3, 1);
%--------------------------------

figure; imshowpair(im, refim);  title('noregistration'); 
figure; imshowpair(newim3, refim); title('ssd');
figure;imshowpair(newim2, refim); title('RC');
figure;imshowpair(newim1, refim);  title('DTV');

% figure; imshowpair(warpAffine2(im_c,res3.X),refim,'blend'); title('ssd');


%  
% % [res2, newim2]=mirt2D_register(refim_dis,im, main, optim);
% I_RMSE3 = RMSE(newim3,refim);
% I_RMSE2 = RMSE(newim2,refim);
% % T_RMSE2 = RMSE(res2.X,T);
% % 
% I_RMSE1 = RMSE(newim1,refim);
% T_RMSE1 = RMSE(res1.X,T);
unregistered = imread('westconcordaerial.png');
figure, imshow(unregistered); 
% text(size(unregistered,2),size(unregistered,1)+15, ...
%     'Image courtesy of mPower3/Emerge', ...
%     'FontSize',7,'HorizontalAlignment','right');
ortho = imread('westconcordorthophoto.png');
figure, imshow(ortho)
% text(size(ortho,2),size(ortho,1)+15, ...
%     'Image courtesy of Massachusetts Executive Office of Environmental Affairs', ...
%     'FontSize',7,'HorizontalAlignment','right');
load westconcordpoints;
input_points = movingPoints;
base_points = fixedPoints;
t_concord = cp2tform(input_points,base_points,'projective');
info = imfinfo('westconcordorthophoto.png');
registered = imtransform(unregistered,t_concord,...
                         'XData',[1 info.Width], 'YData',[1 info.Height]);
figure, imshowpair(registered,ortho,'blend');title('MATLAB');


regsterRC = warpAffine2(im_c,res2);
regsterRC = regsterRC(1:size(ortho,1),1:size(ortho,2),:);
regsterDTV = warpAffine2(im_c,res1);
regsterDTV = regsterDTV(1:size(ortho,1),1:size(ortho,2),:);
figure; imshowpair(regsterRC,refim(1:size(ortho,1),1:size(ortho,2),:),'blend'); 
figure; imshowpair(regsterDTV,refim(1:size(ortho,1),1:size(ortho,2),:),'blend'); 


y = 100:140; x = 240:280;
figure, imshow(unregistered(240:280,77:117,:));
figure, imshowpair(registered(x,y,:),ortho(x,y,:),'blend'); title('matlab');
figure; imshowpair(regsterRC(x,y,:),refim(x,y,:),'blend');  title('RC');
figure; imshowpair(regsterDTV(x,y,:),refim(x,y,:),'blend'); title('DTV');

y = 280:320; x = 150:190;
figure, imshow(unregistered(120:160,258:298,:));
% figure, imshow(refim(x,y,:));
figure, imshowpair(registered(x,y,:),ortho(x,y,:),'blend');title('matlab');
figure; imshowpair(regsterRC(x,y,:),refim(x,y,:),'blend'); title('RC');
figure; imshowpair(regsterDTV(x,y,:),refim(x,y,:),'blend'); title('DTV');