function [raslpara] = ralspara(Transformation)



%% define parameters

% dispaly flag
raslpara.DISPLAY = 0 ;

% save flag
raslpara.saveStart = 1 ;
raslpara.saveEnd = 1 ;
raslpara.saveIntermedia = 0 ;


% for face images
raslpara.canonicalImageSize = [ 120 90];
raslpara.canonicalEyeCoords = [ 16  76 ; ...
                                50 50];
                            
% parametric tranformation model
raslpara.transformType = Transformation; 
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