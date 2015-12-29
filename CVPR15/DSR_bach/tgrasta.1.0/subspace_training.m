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

function[ U_hat, status, Tau_out ] = subspace_training( U_hat, status, FileNames, numtrainingImages, Tau_in, OPTIONS, AlignPara)
%% Using the batch model to train a initial subspace for fully online or simple online model 
QUIET               = OPTIONS.QUIET;
CONVERGE_LEVEL      = OPTIONS.CONVERGE_LEVEL;
subSampling         = OPTIONS.SUBSAMPLING;
max_K               = OPTIONS.NUM_SUBSPACE;
max_outer_loop      = OPTIONS.MAX_OUTER_LOOP;
imgSize             = AlignPara.canonicalImageSize;


%load the images into matrix form
trainingImgs       = cell(numtrainingImages,1);
sigma0 = 2/5 ;  % used for gaussian smoothing
sigmaS = 1 ;    % used for gaussian smoothing
for i=1:numtrainingImages,
  
   currentImage = double(imread(FileNames{i}));
    
    %only use the gray color of each image
    if size(currentImage,3) > 1,   
        currentImage = currentImage(:,:,2);  
    end
    
    currentImage = gamma_decompress(currentImage, 'linear');
    currentImagePyramid = gauss_pyramid( currentImage,1,...   % only consider one level pyramid so far
        sqrt(det(Tau_in{i}(1:2,1:2)))*sigma0, sigmaS );

    trainingImgs{i} = currentImagePyramid{1}; % only consider mono-scale so far
end


%% training subspace

%preparing parameters
Tau_out             = cell(numtrainingImages,1);
DELTA_norms         = nan(numtrainingImages,1);
Jacobians           = cell(numtrainingImages,1);
LocallyalignedImgs  = cell(numtrainingImages,1);
D_aligned           = cell(numtrainingImages,1);
outliers            = cell(numtrainingImages,1);
UwImgs              = cell(numtrainingImages,1);

DIM = imgSize(1)*imgSize(2);

t0 = tic;
for k = 1:max_K, % outerloop for delta_tau update
    fprintf('Delta_tau update %d\n', k);
    
    % prepare the initial image matrix and Jacobian
    for i=1 : numtrainingImages,
        I0_smooth = trainingImgs{i};
        T_in = Tau_in{i};
        
        I0x = imfilter( I0_smooth, (-fspecial('sobel')') / 8 );
        I0y = imfilter( I0_smooth,  -fspecial('sobel')   / 8 );
        
        Tfm = fliptform(maketform('projective',T_in'));
        
        I   = vec(imtransform(I0_smooth, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
        Iu  = vec(imtransform(I0x,Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
        Iv  = vec(imtransform(I0y,Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
        
        y   = I; %vec(I);
        
        Iu = (1/norm(y))*Iu - ( (y'*Iu)/(norm(y))^3 )*y ;
        Iv = (1/norm(y))*Iv - ( (y'*Iv)/(norm(y))^3 )*y ;
        
        % transformation matrix to parameters
        xi = projective_matrix_to_parameters(AlignPara.transformType,T_in) ;
        
        % Compute Jacobian
        Jacobians{i} = image_Jaco(Iu, Iv, imgSize, AlignPara.transformType, xi);
        LocallyalignedImgs{i} = y/norm(y);
    end             

    %plot the inital images 
    if  ~QUIET && k == 1,
        plotTitle = ['initial images'];
        tgrasta_plot_bigpic(LocallyalignedImgs,imgSize,plotTitle,true);
    end

    if k>1,
        U_hat{k} = U_hat{k-1}; status{k} = status{k-1};
        status{k}.level = 2;
    end
    
    tic_inner = tic;
    for inneriter = 1:max_outer_loop,
        image_order = randperm(numtrainingImages);
        
        if status{k}.level >= CONVERGE_LEVEL 
            break;
        end
               
        for i = 1 : numtrainingImages,
            I1 = LocallyalignedImgs{image_order(i)};
            J = Jacobians{image_order(i)};
            tau = Tau_in{image_order(i)};
            
            % subsampling the pre-aligned image, here we subsample the full information           
            p = randperm(DIM);
            idx = p(1: ceil(subSampling * DIM));
            
            
            [ w, e, dual, tau_new, Inew, delta_norm ] = tgrasta_delta_tau_pursuit( U_hat{k}, I1, tau, J, status{k}.OPTS ,AlignPara, idx);
                        
            [ U_hat{k}, status{k}] = tgrasta_update( U_hat{k}, status{k}, w, e, dual, Inew, OPTIONS, idx);
            
            
            UwImgs{image_order(i)} = U_hat{k} * w;
            D_aligned{image_order(i)} = Inew;
            e0 = zeros(DIM,1);
            e0(idx) = e;
            outliers{image_order(i)} = e0;
   
            Tau_out{image_order(i)} = tau_new;
            DELTA_norms(image_order(i)) = delta_norm;
            
        end
       
    end
    
    toc_inner = toc(tic_inner);
    Tau_in = Tau_out;
     
    
    fprintf('Delta_tau %d\n', k);
    fprintf('Outer loop [max delta_tau %.2e, median delta_tau %.2e]\n',...
        max(DELTA_norms), median(DELTA_norms));
    
    % Just training U_hat{1} - U_hat{k}
    if k < max_K && max(DELTA_norms) <= OPTIONS.OUTER_TOL,
        
        for i = (k + 1) : max_K,
            U_hat{i} = U_hat{k};
            status{i} = status{k};
        end
        
        break;    
    end 
    
end
% plot the aligned images
if  ~QUIET ,
        plotTitle = ['aligned images'];
        tgrasta_plot_bigpic(D_aligned,imgSize,plotTitle,true);
end
    
t_training = toc(t0);
fprintf('Training %d images costs %.2f seconds\n', numtrainingImages, t_training);

end

