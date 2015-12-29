% Yigang Peng, Arvind Ganesh, November 2009.
% Questions? abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: RASL: Robust Alignment by Sparse and Low-rank Decomposition for Linearly Correlated Images
%            Yigang Peng, Arvind Ganesh, John Wright, Wenli Xu, and Yi Ma. Proc. of CVPR, 2010.
%

% robust batch image alignment example

% clear
clc ;
clear all;
close all ;

% addpath
addpath RASL_toolbox ;
% addpath data ;
addpath results ;

%% define images' path

subject = 'NUTS';
% FRUITS MOVI NUTS TOY

path = ['datasets/' subject];
Files = dir(strcat(path,'/*.ppm'));
% Files = dir(strcat(path,'\*.pgm'));
LengthFiles = length(Files);
m = 128;
for i=1:LengthFiles
    %         filename = [path  '\' lower(subject)  num2str(i)  '.pgm'];
    filename = [path  '/' lower(subject)  num2str(i)  '.ppm'];
    im = double(imread(filename));
    
    im = mean(im,3);
    im = imresize(im,[m m]);
    %     im = 255*im/max(im(:));
    I0(:,:,i) = im;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = 20;
currentPath = cd;

% im = double(imread('lenargb.jpg'));
% im = mean(im,3);
% im = imresize(im,[m,m]);
% I = repmat(im,[1,1,T]);
[m,n,numImages] = size(I0);

for t=1:numImages
    Itrans = [ 1  0  0
        0  1  0];
    transformations{t} = Itrans;
    
end

tform = [ 1  0  0
    0  1  0
    0 0  1 ];




% output path
destRoot = fullfile(currentPath,'results') ;
destDir = fullfile(destRoot,'lena') ;
if ~exist(destDir,'dir')
    mkdir(destRoot,'lena') ;
end

%% define parameters

% dispaly flag
raslpara.DISPLAY = 0 ;

% save flag
raslpara.saveStart = 1 ;
raslpara.saveEnd = 1 ;
raslpara.saveIntermedia = 0 ;


% for face images
raslpara.canonicalImageSize = [ m m ];
% raslpara.canonicalEyeCoords = [ 5  55 ; ...
%                                 32 32  ];

% parametric tranformation model
raslpara.transformType = 'TRANSLATION';
% one of 'TRANSLATION', 'SIMILARITY', 'AFFINE','HOMOGRAPHY'

raslpara.numScales = 2 ; % if numScales > 1, we use multiscales

% main loop
raslpara.stoppingDelta = .01; % stopping condition of main loop
raslpara.maxIter = 25; % maximum iteration number of main loops

% inner loop
raslpara.inner_tol = 1e-6 ;
raslpara.inner_maxIter = 1000 ;
raslpara.continuationFlag = 1 ;
raslpara.mu = 1e-3 ;
raslpara.lambdac = 1 ; % lambda = lambdac/sqrt(m)


raslpara.alpha = 0.001;
raslpara.beta = 0.0001;
%%
AlignPara.canonicalImageSize = [ m n  ];
AlignPara.transformType = 'TRANSLATION'; % parametric tranformation model,one of 'TRANSLATION', 'EUCLIDEAN', 'SIMILARITY', 'AFFINE','HOMOGRAPHY'
% parameters for subspace update
OPTIONS.NUM_SUBSPACE        = 5; % The number of union of subspaces to approximate the nonlinear transform
OPTIONS.RANK                = 10;
OPTIONS.SUBSAMPLING         = 1.0; % currently subsamping full information
OPTIONS.CONSTANT_STEP       = 0;   % nonzero: constant step

% parameters for sparse_residual_pursuit
OPTIONS.RHO                 = 2;
OPTIONS.ITER_MIN            = 5;
OPTIONS.ITER_MAX            = 20;
OPTIONS.INNER_TOL           = 1e-6;

% batch model
OPTIONS.CONVERGE_LEVEL      = 9;
OPTIONS.OUTER_TOL           = 1e-2;
OPTIONS.MAX_OUTER_LOOP      = 15; % If we want to train the initial subspace
% as clean as possible, it should be large enough!

% simple online model
OPTIONS.MAX_ALIGN_LOOP      = 30;

% simple online model/online model
OPTIONS.stoppingDelta       = 1e-3;

OPTIONS.QUIET               = 1; % 0 display, nonzero no display
OPTIONS.SHOW_INTERVAL       = 1;
%%
for ii=1:1       
    for t=1:numImages
        T(:,:,t) = tform;
        T(3,1,t)  =  round(-5 + 10*rand(1));
        T(3,2,t) =   round(-5 + 10*rand(1));
        
        form_translate = maketform('affine',T(:,:,t));
        I(:,:,t)= imtransform(I0(:,:,t),form_translate, 'XData', [1 (n)],...
            'YData', [1 m]);
        
        Ig(:,:,t)= imtransform(grad(I0(:,:,t)),form_translate, 'XData', [1 (n)],...
            'YData', [1 m]);
    end
    
    AlignPara.init_trans = transformations;
    % profile on;
    cropsize = 0;
    I = I(cropsize+1:end-cropsize,cropsize+1:end-cropsize,:);
    Ig = Ig(cropsize+1:end-cropsize,cropsize+1:end-cropsize,:);
    [m,n,numImages] = size(I);
    ID = reshape(I,[m*n,numImages]);
    figure;showallimages(reshape(I,[m*n,numImages]),[m n]); title('inputs');
    
    t0 = tic;
    [D, Do, A, E, xi, numIterOuter, numIterInner ] = dsr_main_alm(Ig, I, transformations, numImages, raslpara, destDir);
    t1 =   toc;
    
    t0 = tic;
    [D2, Do2, A2, E2, xi2, numIterOuter, numIterInner] = rasl_main2(I, transformations,numImages, raslpara, destDir);
    t2 =   toc;
    
    t0 = tic;
    [UwImgs, alignedImgs, unAlignedImgs, outlierImgs, xi3 ] = tgrasta_batch_training2( I, OPTIONS, AlignPara );
    t3 =   toc;
    
    for jj =1:numImages
        Do3(:,jj) = alignedImgs{jj};
        E3(:,jj) = outlierImgs{jj};
        A3(:,jj)  = UwImgs{jj};
    end
    
    TrueT = -T;
    

%     figure;showallimages(reshape(I0,[m*n,numImages]),[m n]); title('ori');
%     figure;showallimages(D2,raslpara.canonicalImageSize,3); title('inputs');
%     figure;showallimages(Do,raslpara.canonicalImageSize,3); title('DSR');
%     figure;showallimages(Do2,raslpara.canonicalImageSize,3); title('RASL');
%     figure;showallimages(Do3,raslpara.canonicalImageSize,3); title('t-GRASTA');
%     
%     figure;showallimages(E,raslpara.canonicalImageSize,3)
%     figure;showallimages(E2,raslpara.canonicalImageSize,3);
%     figure;showallimages(E3,raslpara.canonicalImageSize,3);
%     
%     figure;showallimages(A,raslpara.canonicalImageSize,3)
%     figure;showallimages(A2,raslpara.canonicalImageSize,3);
%     figure;showallimages(A3,raslpara.canonicalImageSize,3);
    
    
    figure;
    subplot(2,2,1); showallimages(mean(D,2),raslpara.canonicalImageSize,5); title('inputs');
    subplot(2,2,2); showallimages(mean(Do,2),raslpara.canonicalImageSize,5); title('DSR');
    subplot(2,2,3); showallimages(mean(Do2,2),raslpara.canonicalImageSize,5); title('RASL');
    subplot(2,2,4); showallimages(mean(Do3,2),raslpara.canonicalImageSize,5); title('t-GRASTA');
    
    figure;
    subplot(2,2,2); showallimages(mean(A,2),raslpara.canonicalImageSize,5); title('DSR');
    subplot(2,2,3); showallimages(mean(A2,2),raslpara.canonicalImageSize,5); title('RASL');
    subplot(2,2,4); showallimages(mean(A3,2),raslpara.canonicalImageSize,5); title('t-GRASTA');
    
    figure;
    subplot(2,2,2); showallimages((E),raslpara.canonicalImageSize,5); title('DSR');
    subplot(2,2,3); showallimages((E2),raslpara.canonicalImageSize,5); title('RASL');
    subplot(2,2,4); showallimages((E3),raslpara.canonicalImageSize,5); title('t-GRASTA');
    
    %% plot the results
end



