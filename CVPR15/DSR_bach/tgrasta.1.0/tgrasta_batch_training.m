%% t-GRASTA (TRANSFORMED GRASTA)
% Copyright (C) 2013-2014 by Jun He, Dejiao Zhang, and Laura Balzano
%  
%  This file is part of the t-GRASTA library.
%  It is provided without any warranty of fitness
%  for any purpose. You can redistribute this file
%  and/or modify it under the terms of the GNU
%  Lesser General Public License (LGPL) as published
%  by the Free Software Foundation, either version 3
%  of the License or (at your option) any later version.
%  (see http://www.opensource.org/licenses for more info

function [  UwImgs, alignedImgs, unalignedImgs, outliers ] = tgrasta_batch_training( rootPath, trainingDatabaseName, OPTIONS, AlignPara )
%% batch model of t-GRASTA
% If you want to get vary highly aligned frames, you can increase the 
% number of subspace and decrease the value of the stopping creterion-- 
% OPTIONS.OUTER_TOL.            

QUIET               = OPTIONS.QUIET;
CONVERGE_LEVEL      = OPTIONS.CONVERGE_LEVEL;
SHOW_INTERVAL       = OPTIONS.SHOW_INTERVAL;
subSampling         = OPTIONS.SUBSAMPLING;
max_outer_loop      = OPTIONS.MAX_OUTER_LOOP;
max_K               = OPTIONS.NUM_SUBSPACE;

if ~QUIET,
    h_Initimg = subplot(2,2,1);set(gca,'nextplot','replacechildren');title('Locally aligned--(k-1)');
    h_aligimg = subplot(2,2,2);set(gca,'nextplot','replacechildren');title('Locally aligned--(k)');
    h_training_uw = subplot(2,2,3);set(gca,'nextplot','replacechildren');title('I=U*w');
    h_outlier = subplot(2,2,4);set(gca,'nextplot','replacechildren');title('Outlier');
end

baseCoords          = AlignPara.canonicalCoords;
initCoords          = AlignPara.originalCoords;
imgSize             = AlignPara.canonicalImageSize;
DIM                 = imgSize(1)*imgSize(2);

U_hat = cell(max_K,1);    
status.kstatus = cell(max_K,1);
for k=1:max_K,
    status.kstatus{k} = tgrasta_init_status(OPTIONS,DIM);
    U_hat{k} = orth(randn(DIM,OPTIONS.RANK));
end  
init_tau = TwoPointSimilarity( baseCoords, initCoords);
init_tau = [init_tau; 0 0 1] ;


%load the images into matrix form
userDirectoryContents = list_image_files(fullfile(rootPath, trainingDatabaseName));
if isempty(userDirectoryContents)
    error(['No image files were found! Check your paths; there should be images in ' fullfile(rootPath, trainingDatabaseName)]);
end
numtrainingImages = length(userDirectoryContents);

trainingImgs      = cell(numtrainingImages,1);
Tau_in            = cell(numtrainingImages,1);
sigma0 = 2/5 ;  % used for gaussian smoothing
sigmaS = 1 ;    % used for gaussian smoothing
for i=1:numtrainingImages,
  
    frameName = userDirectoryContents{i};
    filename = fullfile(rootPath, trainingDatabaseName, frameName);
    currentImage = double(imread(filename));
      
    Tau_in{i} = init_tau;

    if size(currentImage,3) > 1,   currentImage = currentImage(:,:,2);  end
    
    currentImage = gamma_decompress(currentImage, 'linear');
    
    currentImagePyramid = gauss_pyramid( currentImage,1,...   % only consider one level pyramid so far
    sqrt(det(init_tau(1:2,1:2)))*sigma0, sigmaS );

    trainingImgs{i} = currentImagePyramid{1}; % only consider mono-scale so far
       
end


%% training subspace

Tau_out = cell(numtrainingImages,1);
Jacobians = cell(numtrainingImages,1);
unalignedImgs = cell(numtrainingImages,1);
LocallyalignedImgs = cell(numtrainingImages,1);
alignedImgs = cell(numtrainingImages,1);
outliers    = cell(numtrainingImages,1);
UwImgs  = cell(numtrainingImages,1);
DELTA_norms   = nan(numtrainingImages,1);

% outerloop for delta_tau update
t0 = tic;
for k = 1:max_K, 
%     fprintf('Delta_tau update %d\n', k);
    tic_0=tic;
    %compute transformed images and corresponding Jacobian for delta_tau update
    for i=1 : numtrainingImages,
        I0_smooth = trainingImgs{i};
        Tk_init = Tau_in{i};
        
        I0x = imfilter( I0_smooth, (-fspecial('sobel')') / 8 );
        I0y = imfilter( I0_smooth,  -fspecial('sobel')   / 8 );
        
        Tfm = fliptform(maketform('projective',Tk_init'));
        
        I   = vec(imtransform(I0_smooth, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
        Iu  = vec(imtransform(I0x,Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
        Iv  = vec(imtransform(I0y,Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
        
        y   = I; 
        
        Iu = (1/norm(y))*Iu - ( (y'*Iu)/(norm(y))^3 )*y ;
        Iv = (1/norm(y))*Iv - ( (y'*Iv)/(norm(y))^3 )*y ;
        
        % transformation matrix to parameters
        xi = projective_matrix_to_parameters(AlignPara.transformType,Tk_init) ;
        
        % Compute Jacobian
        Jacobians{i} = image_Jaco(Iu, Iv, imgSize, AlignPara.transformType, xi);
        
        if k == 1,
           unalignedImgs{i} = y/norm(y);
        end
        LocallyalignedImgs{i} = y/norm(y);
    end
    toc_prepare = toc(tic_0);
%     fprintf('prepare aligned images: %.2f sec\n',toc_prepare); 
    
    % initializing the current U_hat{k} with the former value, for batch
    % model we don't necessarily need this step
    if k>1,
        U_hat{k} = U_hat{k-1}; status.kstatus{k} = status.kstatus{k-1};
        status.kstatus{k}.level = 2;
    end
    
    tic_inner = tic;
    for inneriter = 1:max_outer_loop,
        p = randperm(numtrainingImages);
        image_order = p(1:numtrainingImages);
        
        if status.kstatus{k}.level >= CONVERGE_LEVEL , % judge whether the subspace has been comparatively stable/good
            break;
        end
               
        for i = 1 : numtrainingImages,
            
            I1 = LocallyalignedImgs{image_order(i)};
            J = Jacobians{image_order(i)};
            tau = Tau_in{image_order(i)};
            
            % subsampling the pre-aligned image, here we subsample the full information            
            p = randperm(DIM);
            idx = p(1: ceil(subSampling * DIM));
            
            [ w, e, dual, tau_new, Inew, delta_norm ] = tgrasta_delta_tau_pursuit( U_hat{k}, I1, tau, J, status.kstatus{k}.OPTS ,AlignPara, idx);
                        
            [ U_hat{k}, status.kstatus{k}] = tgrasta_update( U_hat{k}, status.kstatus{k}, w, e, dual, Inew, OPTIONS, idx);
            
            Tau_out{image_order(i)} = tau_new;
            
            UwImgs{image_order(i)} = U_hat{k} * w;
            
            alignedImgs{image_order(i)} = Inew;
            
%             e0 = zeros(DIM,1);
%             e0(idx) = e;
            e0 = Inew - U_hat{k} * w;
            outliers{image_order(i)} = e0;
            
            DELTA_norms(image_order(i)) = delta_norm;
            
            if mod(i,SHOW_INTERVAL)==0 && ~QUIET,
                % 1. initial unaligned image
                init_img = reshape(I1,imgSize(1),imgSize(2));
                axes(h_Initimg); imagesc(init_img);colormap gray;axis off;axis ij ;
                
                % 2. current aligned image
                alig_img = reshape( Inew ,imgSize(1),imgSize(2));
                axes(h_aligimg); imagesc(alig_img);colormap gray;axis off;axis ij ;
                
                % 3. lowrank reconstruction
                Uw_img = reshape(U_hat{k} * w, imgSize(1),imgSize(2));
                axes(h_training_uw); imagesc(Uw_img);colormap gray;axis off;axis ij ;
                
                % 4. outliers
                e_img = zeros(DIM,1); e_img(idx) = e;
                e_img = reshape(e_img,imgSize(1),imgSize(2));
                axes(h_outlier); imagesc(e_img);colormap gray;axis off;axis ij ;
                
            end           
        end   
    end
    % update tau
    Tau_in = Tau_out;
    
    toc_inner = toc(tic_inner);
%     fprintf('Delta_tau %d: %.2f seconds\n', k, toc_inner);
    fprintf('Delta_tau %d:\n', k);

    fprintf('Outer loop [max delta_tau %.2e, median delta_tau %.2e ,tol %.2e]\n',...
        max(DELTA_norms), median(DELTA_norms),OPTIONS.OUTER_TOL);
    
    if k == max_K || max(DELTA_norms) < OPTIONS.OUTER_TOL        
        t_training = toc(t0);
        fprintf('Using %d subspaces to training %d images costs %.2f seconds\n', k, numtrainingImages, t_training);
        break;
    end
end

end

