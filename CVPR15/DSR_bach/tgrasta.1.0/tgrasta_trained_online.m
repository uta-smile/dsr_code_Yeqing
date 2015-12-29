% t-GRASTA (TRANSFORMED GRASTA)
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

function [  UwImgs, alignedImgs, unAlignedImgs, outlierImgs ] = tgrasta_trained_online( rootPath, trainingDatabaseName, OPTIONS, AlignPara )
% Trained online model of t-GRASTA
% We first randomly select a small part of the tatal images and use the
% bathch model of t-GRASTA to train a initial subspace.
% And then we use this trained online mode to align the rest images without updating the 
% subspace.
%
% Note: we usually randomly select 10% - 30% (depended on different
% datasets and requirements) of the total images to train
% the initial subspace in order to get a "well-trained " subspace
%
if ~OPTIONS.QUIET,
    h_Initimg = subplot(2,2,1);set(gca,'nextplot','replacechildren');title('Unaligned');
    h_aligimg = subplot(2,2,2);set(gca,'nextplot','replacechildren');title('Aligned');
    h_training_uw = subplot(2,2,3);set(gca,'nextplot','replacechildren');title('I=U*w');
    h_outlier = subplot(2,2,4);set(gca,'nextplot','replacechildren');title('Outlier');
end

SHOW_INTERVAL       = OPTIONS.SHOW_INTERVAL;
max_K               = OPTIONS.NUM_SUBSPACE;
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


% getting images 
userDirectoryContents = list_image_files(fullfile(rootPath, trainingDatabaseName));
if isempty(userDirectoryContents)
    error(['No image files were found! Check your paths; there should be images in ' fullfile(rootPath, trainingDatabaseName)]);
end
numImages           = length(userDirectoryContents);


%% obtaining a well-trained subspace 

% randomly selecting a small part of total images 
numtrainingImages = min(150, 0.20 * numImages); 
p = randperm(numImages);
image_sampled = p(1:numtrainingImages);

Transformation = cell(numImages,1);
transformation = cell(numtrainingImages,1);
FileNames     = cell(numtrainingImages,1);
for i=1:numImages,
    
   Transformation{i} = init_tau;
   
   if i <= numtrainingImages,     
       transformation{i} = init_tau;
       frameName = userDirectoryContents{image_sampled(i)};
       FileNames{i}= fullfile(rootPath, trainingDatabaseName, frameName);    
   end
   
end

[ U_hat, status.kstatus, transformation] = subspace_training( U_hat, status.kstatus, FileNames, numtrainingImages, transformation , OPTIONS, AlignPara);

% storing the transformation obtained by the batch model to serve as the
% inital value for the following simple alignment process
for i = 1 : numtrainingImages,
    Transformation{image_sampled(i)} = transformation{i};
end

    
%% aligning images without subspace update
UwImgs              = cell(numImages,1);
alignedImgs         = cell(numImages,1);
unAlignedImgs       = cell(numImages,1);
outlierImgs         = cell(numImages,1);
U0                  = U_hat{max_K};% well-trained subspace

t0 = tic;
for fileIndex = 1:numImages,
        
        frameName = userDirectoryContents{fileIndex};
        frameFileName = fullfile(rootPath, trainingDatabaseName, frameName);
        Tau = Transformation{i};
         fprintf('frame %d :\n', fileIndex);
         
         [ w, Inew, Iinit, e] = aligning_image( U0, frameFileName, Tau, OPTIONS, AlignPara);
    
         UwImgs{fileIndex} = U0 * w;
         alignedImgs{fileIndex} = Inew;
         unAlignedImgs{fileIndex} = Iinit;
         outlierImgs{fileIndex} = Inew - U0 * w;
%          outlierImgs{fileIndex} = e;

    if mod(fileIndex,SHOW_INTERVAL)==0 && ~OPTIONS.QUIET,
        % 1. initial unaligned image
        init_img = reshape(Iinit,imgSize(1),imgSize(2));
        axes(h_Initimg); imagesc(init_img);colormap gray;axis off;axis ij ;

        % 2. current aligned image
        alig_img = reshape( Inew ,imgSize(1),imgSize(2));
        axes(h_aligimg); imagesc(alig_img);colormap gray;axis off;axis ij ;

        % 3. lowrank reconstruction
        Uw_img = reshape(U0 * w, imgSize(1),imgSize(2));
        axes(h_training_uw); imagesc(Uw_img);colormap gray;axis off;axis ij ;

        % 4. outliers
        e = Inew - U0 * w;
        e_img = reshape(e,imgSize(1),imgSize(2));
        axes(h_outlier); imagesc(e_img);colormap gray;axis off;axis ij ;

    end
end
        
t_online = toc(t0);
fprintf('Online aligning: %.2f seconds\n', t_online);

end


