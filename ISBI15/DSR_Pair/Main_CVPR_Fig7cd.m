% MIRT2D_EXAMPLE2: Non-rigid 2D registration example 2 with Residual Complexity (RC)
% similarity measure, where one of the images is corrupted by nonstationary
% intensity distortion. Try using other similarity measures, e.g. SSD, CC
% or MI to see that they won't work with nonstationary instensity distortions

clear all; close all; clc;
refim = imread('lenargb.jpg');

[m,n] = size(refim);
refim = double(refim)/255;

refim = imresize(refim,[128,128]);
% I = zeros(150,150,3);
% I(12:end-11,12:end-11,:) = refim;
im = refim;
Tim = [1.00 0.00 0
    -0.00 1.00 0];

im = rgb2gray(im);




% Main settings
main.similarity='dtv';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI
main.type = 'AFFINE'; %translation
main.TV = TVOP;
main.subdivide=3;       % use 3 hierarchical levels
main.single=0;          % show mesh transformation at every iteration
main.alpha = 0.05;

% Optimization settings
optim.maxsteps = 100;   % maximum number of iterations at each hierarchical level
optim.fundif = 1e-8;    % tolerance (stopping criterion)
optim.gamma = 1;       % initial optimization step size
optim.anneal=0.8;       % annealing rate on the optimization step

main2 = main;
main2.similarity = 'rc';
% main.similarity='dtv';

main3= main;
main3.similarity = 'ssd';

% main4 = main;
% main4.similarity = 'djtv';

% im = im_int(:,:,1);
% refim_dis = add_distor(refim,2,0);
% refim_dis = add_distor(refim_dis,2,1);
% refim_dis = refim;


display = 0;

% K = 1:2;Iter=10;
K = 1:6;Iter=2; % it can be changed, 50 runs used in the paper
for k=K
    for i = 1:Iter
        T = Tim + 0.02*rand(size(Tim));
        T(1:2,3) = T(1:2,3) + 2*rand(2,1);
        T = [T;[0,0,1]];
        refim = warpAffine2(im,T);
        
        
        
%         refim_dis = add_distor(refim,1,0);
        refim_dis = add_distor(refim,k);
        
        if display
            figure;
            subplot(2,3,1);imshow(refim_dis); title('Reference (fixed) image');
            subplot(2,3,2);imshow(im); title('Source image');
        end
        
        t0= cputime();
        [res1, newim1]=mirt2D_register_rigid(refim_dis,im, main, optim);
        t1(k,i)= cputime()-t0;
        I_RMSE1(k,i) = RMSE(newim1,refim);
        T_RMSE1(k,i) = RMSE(res1,T);
        
        t0= cputime();
        [res2, newim2]=mirt2D_register_rigid(refim_dis,im, main2, optim);
        t2(k,i)= cputime()-t0;
        I_RMSE2(k,i) = RMSE(newim2,refim);
        T_RMSE2(k,i) = RMSE(res2,T);
        
        t0= cputime();
        [res3, newim3]=mirt2D_register_rigid(refim_dis,im, main3, optim);
        t3(k,i)= cputime()-t0;
        I_RMSE3(k,i) = RMSE(newim3,refim);
        T_RMSE3(k,i) = RMSE(res3,T);
        
%         t0= cputime();
%         [res4, newim4]=mirt2D_register_rigid(refim_dis,im, main4, optim);
%         t4(k,i)= cputime()-t0;
%         I_RMSE4(k,i) = RMSE(newim4,refim);
%         T_RMSE4(k,i) = RMSE(res4.X,T);
        
        if display
            figure,subplot(2,3,1); imshowpair(im, refim_dis);  title('UnRegistered');
            subplot(2,3,2); imshowpair(newim1, refim_dis); title('DTV');
            subplot(2,3,3); imshowpair(newim2, refim_dis); title('RC');
            subplot(2,3,4); imshowpair(newim3, refim_dis); title('SSD');
%             subplot(2,3,5); imshowpair(newim4, refim_dis); title('DJTV');
        end
    end
end


M_RMSE1 = mean(I_RMSE1,2); S_RMSE1 = std(I_RMSE1,1,2);
M_RMSE2 = mean(I_RMSE2,2); S_RMSE2 = std(I_RMSE2,1,2);
M_RMSE3 = mean(I_RMSE3,2); S_RMSE3 = std(I_RMSE3,1,2);
% M_RMSE4 = mean(I_RMSE4,2); S_RMSE4 = std(I_RMSE4,1,2);

MT_RMSE1 = mean(T_RMSE1,2); ST_RMSE1 = std(T_RMSE1,1,2);
MT_RMSE2 = mean(T_RMSE2,2); ST_RMSE2 = std(T_RMSE2,1,2);
MT_RMSE3 = mean(T_RMSE3,2); ST_RMSE3 = std(T_RMSE3,1,2);
% MT_RMSE4 = mean(T_RMSE4,2); ST_RMSE4 = std(T_RMSE4,1,2);


ls=3; ms=8; ts=20;
figure; hold on;box on;
errorbar(K,M_RMSE3, S_RMSE3, 'g--', 'linewidth',ls, 'markersize', ms);
errorbar(K,M_RMSE2, S_RMSE2, 'b-.', 'linewidth',ls, 'markersize', ms);
errorbar(K,M_RMSE1, S_RMSE1, 'r-', 'linewidth',ls, 'markersize', ms);
% errorbar(K,M_RMSE4, S_RMSE4, 'c-', 'linewidth',ls, 'markersize', ms);
legend('SSD','RC','Proposed','DJTV','Location','SouthEast');
xlabel('K');
ylabel('Intensity RMSE');
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts);
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts);
axis([0 7 -0.02 0.16]);
set(gca,'FontSize',16);

% ls=2; ms=8; ts=20;
figure; hold on;box on;
errorbar(K,MT_RMSE3, ST_RMSE3, 'g--', 'linewidth',ls, 'markersize', ms);
errorbar(K,MT_RMSE2, ST_RMSE2, 'b-.', 'linewidth',ls, 'markersize', ms);
errorbar(K,MT_RMSE1, ST_RMSE1, 'r-', 'linewidth',ls, 'markersize', ms);
% errorbar(K,MT_RMSE4, ST_RMSE4, 'c-', 'linewidth',ls, 'markersize', ms);
legend('SSD','RC','Proposed','DJTV','Location','SouthEast');
xlabel('K');
ylabel('Transformation RMSE');
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts);
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts);
axis([0 7 0 2]);
set(gca,'FontSize',16);




Time1 = mean(t1(:))
Time2 = mean(t2(:))
Time3 = mean(t3(:))
% Time4 = mean(t4(:))
