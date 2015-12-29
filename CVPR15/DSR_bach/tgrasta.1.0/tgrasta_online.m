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
%  (see http://www.opensource.org/licenses for more info)

%
function [ U_hat, status, Inew, Iinit ] = tgrasta_online( U_hat, status, filename, T_in, OPTIONS,AlignPara)
%% fully online model of t-GRASTA for one frame, used in the function tgrasta_training_online
%

MAX_INNER_ITER   = 1;
max_K            = OPTIONS.NUM_SUBSPACE;
subSampling      = OPTIONS.SUBSAMPLING;

imgSize          = AlignPara.canonicalImageSize;
DIM              = imgSize(1)*imgSize(2);

sigma0 = 2/5 ;  % used for gaussian smoothing
sigmaS = 1 ;    % used for gaussian smoothing
currentImage = double(imread(filename)); 
if size(currentImage,3) > 1,  %only use the gray color of each image 
    currentImage = currentImage(:,:,2);  
end 
currentImage = gamma_decompress(currentImage, 'linear');
currentImagePyramid = gauss_pyramid( currentImage,1,...   % only consider one level pyramid so far
    sqrt(det(T_in(1:2,1:2)))*sigma0, sigmaS );

I0 = currentImagePyramid{1}; % only consider mono-scale so far
I0x = imfilter( I0, (-fspecial('sobel')') / 8 );
I0y = imfilter( I0,  -fspecial('sobel')   / 8 );
for k=1:max_K,

    Tfm = fliptform(maketform('projective',T_in'));
    I   = vec(imtransform(I0, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
    Iu  = vec(imtransform(I0x,Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
    Iv  = vec(imtransform(I0y,Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
    
    Iu = (1/norm(I))*Iu - ( (I'*Iu)/(norm(I))^3 )*I ;
    Iv = (1/norm(I))*Iv - ( (I'*Iv)/(norm(I))^3 )*I ;
    
    % transformation matrix to parameters
    xi = projective_matrix_to_parameters(AlignPara.transformType,T_in) ;
    
    % Compute Jacobian
    J = image_Jaco(Iu, Iv, imgSize, AlignPara.transformType, xi);
    I1 = I/norm(I);

    if k==1, Iinit = I1; end
    
    % subsampling the pre-aligned image
    p = randperm(DIM);
    idx = p(1: ceil(subSampling * DIM));
    
   
    converge = false; 
    inner_k = 0;
    while ~converge && inner_k < MAX_INNER_ITER,
        inner_k = inner_k +1;
        
        [ w, e, dual, T_out, Inew, delta_norm ] = tgrasta_delta_tau_pursuit( U_hat{k}, I1, T_in, J, status{k}.OPTS ,AlignPara, idx);
        
        [ U_hat{k}, status{k} ] = tgrasta_update( U_hat{k}, status{k}, w, e, dual, Inew, OPTIONS,idx);
       
        if delta_norm < OPTIONS.stoppingDelta,
            converge = true;
        end
        
    end
    T_in = T_out;
end

end

