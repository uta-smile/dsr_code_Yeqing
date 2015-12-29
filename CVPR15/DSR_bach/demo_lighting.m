
% clear
clc ;
clear all;
close all ;

% addpath
addpath RASL_toolbox ;
addpath data ;
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

currentPath = cd;
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
raslpara.numScales = 2; % if numScales > 1, we use multiscales
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
% parameters for image alignment
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
I(:,:,1)=I0(:,:,1);
Ig(:,:,1)=grad(I0(:,:,1));
%%
for ii=1:1  
    for t=2:numImages
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
    [UwImgs, alignedImgs, unAlignedImgs, outlierImgs, xi3 ] = tgrasta_batch_training2( I, OPTIONS, AlignPara );
     
    [D, Do, A, E, xi, numIterOuter, numIterInner] = rasl_main2(I, transformations,numImages, raslpara, destDir);
   
    [D2, Do2, A2, E2, xi2, numIterOuter, numIterInner ] = dsr_main_alm(Ig, I, transformations, numImages, raslpara, destDir);
          
    TrueT = -T;

    for jj =2:numImages
    
    xi{jj} = xi{jj}  - xi{1};
    xi2{jj} = xi2{jj} - xi2{1};    
    xi3{jj} = xi3{jj} - xi3{1};      
        
    err0(ii,jj) = mean(abs(T(3,1:2,jj)));
    err1(ii,jj) = mean(abs(xi{jj}'-T(3,1:2,jj)));
    err2(ii,jj) = mean(abs(xi2{jj}'-T(3,1:2,jj)));
    err3(ii,jj) = mean(abs(xi3{jj}(1:2,3)'-T(3,1:2,jj)));
    end
     
    %% plot the results
end
% disp('mean.........')
% mean(err0(:))
% mean(err1(:))
% mean(err2(:))
% mean(err3(:))
% disp('max...........')
% max(err0(:))
% max(err1(:))
% max(err2(:))
% max(err3(:))

fprintf('Unregistered: mean error: %f, max error: %f \n',mean(err0(:)),max(err0(:)));
fprintf('tgrasta: mean error: %f, max error: %f \n',mean(err3(:)),max(err3(:)));
fprintf('RASL: mean error: %f, max error: %f \n',mean(err1(:)),max(err1(:)));
fprintf('DSR: mean error: %f, max error: %f \n',mean(err2(:)),max(err2(:)));




