
% clear
clc ;
clear all;
close all ;

% addpath
addpath RASL_toolbox ;
addpath data ;
addpath results ;

% for sub = 1:100
 for sub = 1:1
    %% define images' path
    
    subject = num2str(sub,'%03d');
    
    path = ['P:\data\multiview\' subject '\01\05_1'];
    Files = dir(strcat(path,'\*.png'));
    % Files = dir(strcat(path,'\*.pgm'));
    LengthFiles = length(Files);
    % m = 128;
    for i=1:LengthFiles
        %         filename = [path  '\' lower(subject)  num2str(i)  '.pgm'];
        %     filename = [path  '\' lower(subject)  num2str(i)  '.png'];
        filename = Files(i).name;
        filename = [path '\' filename];
        im = double(imread(filename));
        
        im = mean(im,3);
        im = imresize(im,[240 320]);
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
    %
    
    
    tform = [ 1  0  0
        0  1  0];
    %         0 0  1 ];
    
    
    
    
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
    raslpara.canonicalImageSize = [120 90];
    raslpara.canonicalEyeCoords = [ 16  76 ; ...
        50 50];
    
    % parametric tranformation model
    raslpara.transformType = 'EUCLIDEAN';
    % one of 'TRANSLATION', 'SIMILARITY', 'AFFINE','HOMOGRAPHY'  EUCLIDEAN
    
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
    % AlignPara.originalCoords = [20 60; ...
    %                        49.5 49.5];
    AlignPara.canonicalImageSize = raslpara.canonicalImageSize;
    %     AlignPara.canonicalCoords = raslpara.canonicalEyeCoords;
    
    
    AlignPara.transformType = raslpara.transformType; % parametric tranformation model,one of 'TRANSLATION', 'EUCLIDEAN', 'SIMILARITY', 'AFFINE','HOMOGRAPHY'
    % parameters for subspace update
    OPTIONS.NUM_SUBSPACE        = 10; % The number of union of subspaces to approximate the nonlinear transform
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
    %%% Get training images
    
    % get initial transformation
    % transformationInit = 'SIMILARITY';
    % [fileNames, transformations, numImages] = get_training_images( imagePath, pointPath, userName, raslpara.canonicalEyeCoords, transformationInit) ;
    
    
    
    
    %%
    faceDetector = vision.CascadeObjectDetector;
    bboxes = step(faceDetector, I0(:,:,1)/255);
    m = bboxes(3);n = bboxes(4);
    xstart = bboxes(1);    ystart = bboxes(2);
    
    
    %         [m,n] =
    % T = [];
    %%
    for ii=1:1
        for t=1:numImages
            T(:,:,t) = tform;
            
            T(1,3,t)  = xstart;
            T(2,3,t) = ystart;
            
            T(1,3,t)  =  T(1,3,t) + round(-5 + 10*rand(1));
            T(2,3,t)  =   T(2,3,t)  + round(-5 + 10*rand(1));
            
            
            
            theta = pi*(-5+10*rand(1))/180;
            %             T(1:2,1:2,t)  =  T(1:2,1:2,t) + (-0.05+0.1*rand(2,2));
            
            T(1:2,1:2,t) = [cos(theta) -sin(theta) ; sin(theta)  cos(theta) ];
            
            transformations{t} = T(:,:,t);
            Ig(:,:,t) = grad(I0(:,:,t));
            
        end
        
        AlignPara.init_trans = transformations;
        
        tic;
        [UwImgs, alignedImgs, unAlignedImgs, outlierImgs, XI3 ] = tgrasta_batch_training2( I0, OPTIONS, AlignPara );
        t3 = toc;
        
            tic;
        [D, Do, A, E, xi, numIterOuter, numIterInner ] = dsr_main2(Ig, I0,transformations, numImages, raslpara, destDir);
                t1 = toc;
                
                      tic;
        [Din, Do2, A2, E2, xi2, numIterOuter, numIterInner] = rasl_main2(I0, transformations,numImages, raslpara, destDir);
                   t2 = toc;
        
        for jj =1:numImages
            Do3(:,jj) = alignedImgs{jj};
            E3(:,jj) = outlierImgs{jj};
            A3(:,jj)  = UwImgs{jj};
            
            beta3(:,jj) =  projective_matrix_to_parameters(raslpara.transformType,XI3{jj});
            beta2(:,jj) = xi2{jj};
            beta(:,jj) = xi{jj};
            
        end
        
        std_dsr(ii,:) =  std(beta,1,2);
        std_rasl(ii,:) =  std(beta2,1,2);
        std_tgra(ii,:) =  std(beta3,1,2);
    end
     std_dsr_100(sub,:) =  mean(std_dsr,1);
     std_rasl_100(sub,:) =  mean(std_rasl,1);
     std_tgra_100(sub,:) =  mean(std_tgra,1);
end
ls=2; ms=1; ts=16;
figure;hold on
plot(std_rasl_100(:,1)*57.29,'g-.', 'linewidth',ls, 'markersize', ms);
plot(std_tgra_100(:,1)*57.29,'b--', 'linewidth',ls, 'markersize', ms);
plot(std_dsr_100(:,1)*57.29,'r-', 'linewidth',ls, 'markersize', ms);
xlabel('Subject');ylabel('Rotation STD');box on;
legend('RASL','t-GRASTA','Proposed', 'Location','SouthEast');
textobj = findobj('type', 'text'); set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel'); set(h_xlabel,'FontSize',ts);
h_xlabel = get(gca,'YLabel');set(h_xlabel,'FontSize',ts);

figure;hold on
plot(std_rasl_100(:,2),'g-.', 'linewidth',ls, 'markersize', ms);
plot(std_tgra_100(:,2),'b--', 'linewidth',ls, 'markersize', ms);
plot(std_dsr_100(:,2),'r-', 'linewidth',ls, 'markersize', ms);
xlabel('Subject');ylabel('X-translation STD');box on;
legend('RASL','t-GRASTA','Proposed', 'Location','SouthEast');
textobj = findobj('type', 'text'); set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel'); set(h_xlabel,'FontSize',ts);
h_xlabel = get(gca,'YLabel');set(h_xlabel,'FontSize',ts);

figure;hold on
plot(std_rasl_100(:,3),'g-.', 'linewidth',ls, 'markersize', ms);
plot(std_tgra_100(:,3),'b--', 'linewidth',ls, 'markersize', ms);
plot(std_dsr_100(:,3),'r-', 'linewidth',ls, 'markersize', ms);
xlabel('Subject');ylabel('Y-translation STD');box on;
legend('RASL','t-GRASTA','Proposed', 'Location','SouthEast');
textobj = findobj('type', 'text'); set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel'); set(h_xlabel,'FontSize',ts);
h_xlabel = get(gca,'YLabel');set(h_xlabel,'FontSize',ts);
