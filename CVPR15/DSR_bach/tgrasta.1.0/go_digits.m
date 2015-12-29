%% Need check!
clc ;
clear all;
close all ;

addpath RASL_toolbox_2010 ;
addpath data ;

% input path
currentPath = cd;
imagePath = fullfile(currentPath,'data') ;
TrainingUserName = 'Digits' ;

% parameters for image alignment
AlignPara.originalCoords = [14 24;...
                      18 18 ];

AlignPara.canonicalImageSize = [ 29 29  ];
AlignPara.canonicalCoords = [ 10 20 ; ...
                      14 14  ];
AlignPara.transformType = 'EUCLIDEAN'; 

% parameters for subspace update 
OPTIONS.NUM_SUBSPACE        = 10; %The number of union of subspaces to approximate the nonlinear transform
OPTIONS.RANK                = 15;
OPTIONS.SUBSAMPLING         = 1.0; % currently subsamping full information
OPTIONS.CONSTANT_STEP       = 0;%5*1e-3;   % nonzero: constant step


% parameters for sparse_residual_pursuit
OPTIONS.RHO                 = 2;   
OPTIONS.ITER_MIN            = 5;
OPTIONS.ITER_MAX            = 20;
OPTIONS.INNER_TOL           = 1e-6;

% batch model
OPTIONS.CONVERGE_LEVEL      = 9;   
OPTIONS.OUTER_TOL           = 1e-2; 
OPTIONS.MAX_OUTER_LOOP      = 5; 

% simple online model
OPTIONS.MAX_ALIGN_LOOP      = 30; 

% simple online model/online model
OPTIONS.stoppingDelta       = 1e-3; 

OPTIONS.QUIET               = 1; % 0 display, nonzero no display
OPTIONS.SHOW_INTERVAL       = 1; 


%% Aligning image by t-GRASTA with certain model

MODEL         = 'batch'; % should be one of "online","simple" and "batch"

if strcmp(MODEL,'online'), 
    [ UwImgs, alignedImgs, unAlignedImgs, outlierImgs ] = tgrasta_fully_online( imagePath, TrainingUserName, OPTIONS, AlignPara);  
elseif strcmp(MODEL,'simple'), 
    [ UwImgs, alignedImgs, unAlignedImgs, outlierImgs ] = tgrasta_trained_online( imagePath, TrainingUserName, OPTIONS, AlignPara );
elseif strcmp(MODEL,'batch')
    [ UwImgs, alignedImgs, unAlignedImgs, outlierImgs] = tgrasta_batch_training( imagePath, TrainingUserName, OPTIONS, AlignPara );
end

% %% make seperation video 
% video_name = TrainingUserName;
% vInfo.rows = AlignPara.canonicalImageSize(1);
% vInfo.cols = AlignPara.canonicalImageSize(2);
% make_video( video_name, outlierImgs,UwImgs, alignedImgs, unAlignedImgs,vInfo);

%% Plot results
tgrasta_plot_bigpic(unAlignedImgs, AlignPara.canonicalImageSize , 'Unaligned Images',1, 1);
tgrasta_plot_bigpic(alignedImgs, AlignPara.canonicalImageSize , 'Aligned Images',1, 1);
tgrasta_plot_bigpic(UwImgs, AlignPara.canonicalImageSize , 'Low-rank U*W',1, 1);
tgrasta_plot_bigpic(outlierImgs, AlignPara.canonicalImageSize , 'Outliers e',1, 1);