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

function [ w, Inew, Iinit,e ] = aligning_image( U0, filename, T_in, OPTIONS, AlignPara)
%% Aligning images with a well-trained subspace without subspace update
% Used in the simple online mode.

OPTS.MAX_ITER = OPTIONS.ITER_MAX;
OPTS.TOL      = OPTIONS.INNER_TOL;
OPTS.RHO      = OPTIONS.RHO;
MAX_LOOP      = OPTIONS.MAX_ALIGN_LOOP;
subSampling   = OPTIONS.SUBSAMPLING;
imgSize       = AlignPara.canonicalImageSize;
DIM           = imgSize(1)*imgSize(2);

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

converge = false;
numiter = 0;
 while ~converge && numiter< MAX_LOOP,

    numiter = numiter +1;
    
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
    if numiter == 1, Iinit = I1; end
    
    % subsampling the pre-aligned image
    p = randperm(DIM);
    idx = p(1: ceil(subSampling * DIM));
    
    [ w, e, ~, T_out, Inew, delta_norm ] = tgrasta_delta_tau_pursuit( U0, I1, T_in, J, OPTS ,AlignPara, idx);
    
    if delta_norm < OPTIONS.stoppingDelta ,
        converge = true;
    end
    
    T_in = T_out;
               
 end
     
end



