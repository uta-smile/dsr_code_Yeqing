
% clear
clc ;
clear all;
close all ;
dbstop if error;

% addpath
addpath RASL_toolbox ;
addpath data ;
addpath results ;

%% define images' path
path = '/Users/yeqing/Documents/Data/lfwcrop_grey/faces';
Files = dir(fullfile(path,'*.pgm'));
LengthFiles = length(Files);
m = 64;

fileID = fopen(fullfile('results', 'lfwcrop.txt'),'a');
max_iter = 1;

for num_sample = [50, 80, 100],

subset = randperm(LengthFiles, num_sample);
I0 = zeros(m, m, num_sample);
for i=1:num_sample%1:LengthFiles
    j = subset(i);
    filename = fullfile(path, Files(j).name);
    im = double(imread(filename));    
    im = imresize(im,[m m]);
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
for ii=1:max_iter,
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
     
    tic;
    [D, Do, A, E, xi, numIterOuter, numIterInner] = rasl_main2(I, transformations,numImages, raslpara, destDir);
    t_rasl = toc;
   
    tic;
    [D2, Do2, A2, E2, xi2, numIterOuter, numIterInner ] = dsr_main_alm(Ig, I, transformations, numImages, raslpara, destDir);
    t_dsr = toc;
    
    tic;
    AlignPara.init_trans = transformations;
    [UwImgs, alignedImgs, unAlignedImgs, outlierImgs, xi3 ] = tgrasta_batch_training2( I, OPTIONS, AlignPara );
    t_gra = toc;
          
    TrueT = -T;

    time_rasl(ii) = t_rasl;
    time_dsr(ii) = t_dsr;
    time_gra(ii) = t_gra;
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

fprintf(fileID, '%d/%d samples m=%d\n', num_sample, LengthFiles, m);
fprintf(fileID, 'Unregistered: mean error: %f, max error: %f \n',mean(err0(:)),max(err0(:)));
fprintf(fileID, 'tgrasta: mean error: %f, max error: %f \n',mean(err3(:)),max(err3(:)));
fprintf(fileID, 'RASL: mean error: %f, max error: %f \n',mean(err1(:)),max(err1(:)));
fprintf(fileID, 'DSR: mean error: %f, max error: %f \n',mean(err2(:)),max(err2(:)));

fprintf(fileID, 'tgrasta: mean running time: %f\n',mean(time_gra(:)));
fprintf(fileID, 'RASL: mean running time: %f\n',mean(time_rasl(:)));
fprintf(fileID, 'DSR: mean running time: %f\n',mean(time_dsr(:)));

end

fclose(fileID);

