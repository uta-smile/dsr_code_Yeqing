function [AlignPara,OPTIONS] = tgsrpara(Transformation)



% parameters for image alignment

AlignPara.canonicalImageSize = [ 200 150  ];
AlignPara.canonicalCoords = [31  111; ...
                    105 94];

AlignPara.transformType = Transformation; % parametric tranformation model,one of 'TRANSLATION', 'EUCLIDEAN', 'SIMILARITY', 'AFFINE','HOMOGRAPHY'
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
