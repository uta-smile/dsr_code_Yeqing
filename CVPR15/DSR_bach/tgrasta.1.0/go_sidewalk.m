
clc ;
clear all;
close all ;

addpath RASL_toolbox_2010 ;
addpath data ;

% input path
currentPath = cd;
imagePath = fullfile(currentPath,'data') ;
TrainingUserName = 'Sidewalk' ;

% parameters for image alignment
AlignPara.originalCoords = [90 180;...
                       160 160 ];
AlignPara.canonicalImageSize = [ 50 78  ];
AlignPara.canonicalCoords = [ 17.5 40 ; ...
                      35 35  ];
AlignPara.transformType = 'AFFINE'; 

% parameters for subspace update 
OPTIONS.NUM_SUBSPACE        = 5; %The number of union of subspaces to approximate the nonlinear transform
OPTIONS.RANK                = 5;
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


%% Aligning image by t-GRASTA with certain model

MODEL         = 'online'; % should be one of "online","simple" and "batch"

if strcmp(MODEL,'online'), 
    [ UwImgs, alignedImgs, unAlignedImgs, outlierImgs ] = tgrasta_fully_online( imagePath, TrainingUserName, OPTIONS, AlignPara);  
elseif strcmp(MODEL,'simple'), 
    [ UwImgs, alignedImgs, unAlignedImgs, outlierImgs ] = tgrasta_trained_online( imagePath, TrainingUserName, OPTIONS, AlignPara );
elseif strcmp(MODEL,'batch')
    [ UwImgs, alignedImgs, unAlignedImgs, outlierImgs] = tgrasta_batch_training( imagePath, TrainingUserName, OPTIONS, AlignPara );
end

%% make seperation video 
video_name = TrainingUserName;
vInfo.rows = AlignPara.canonicalImageSize(1);
vInfo.cols = AlignPara.canonicalImageSize(2);
make_video( video_name, outlierImgs,UwImgs, alignedImgs, unAlignedImgs,vInfo);