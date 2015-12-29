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

function [  UwImgs, alignedImgs, unAlignedImgs, outlierImgs ] = tgrasta_fully_online( rootPath, trainingDatabaseName, OPTIONS, AlignPara )
%% Fully online model of t-GRASTA, realted to the function tgrasta_online
%
% We first using batch mode to train a initial subspace by the  first 20
% frames. And then we use it for the following online image alignment 
% process.
%
% Note: for the online model, we can also directly go to the online image
% alignment process without the initial subsapce. However, these is a
% trade-off between the alignment quality of the first dozens of images and
% the whole efficiency
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
numImages =length(userDirectoryContents); % min(length(userDirectoryContents), 600);


%% training the initial subspace using batch model

numtrainingImages = 20;
p = 1:numImages;
image_sampled = p(1:numtrainingImages);

Transformation = cell(numtrainingImages,1);
FileNames     = cell(numtrainingImages,1);
for i=1:numtrainingImages,
    
   Transformation{i} = init_tau; 
   
   frameName = userDirectoryContents{image_sampled(i)};
   FileNames{i}= fullfile(rootPath, trainingDatabaseName, frameName);
   
end

fprintf('Firstly,using %d frames to train the initial subspaces\n',numtrainingImages);

[ U_hat, status.kstatus, ~] = subspace_training( U_hat, status.kstatus, FileNames, numtrainingImages, Transformation , OPTIONS, AlignPara);

%% online image alignment
fprintf('Online aligning %d frames\n',numImages);

UwImgs              = cell(numImages,1);
alignedImgs         = cell(numImages,1);
unAlignedImgs       = cell(numImages,1);
outlierImgs         = cell(numImages,1);

t0 = tic;
for fileIndex = 1:numImages,
        
        frameName = userDirectoryContents{fileIndex};
        frameFileName = fullfile(rootPath, trainingDatabaseName, frameName);
       
         fprintf('frame %d :\n', fileIndex);
        [ U_hat, status.kstatus, Inew, Iinit ] = tgrasta_online( U_hat, status.kstatus, frameFileName, init_tau , OPTIONS, AlignPara);
       
        
        UwImgs{fileIndex} = U_hat{max_K} * status.kstatus{max_K}.w;
        alignedImgs{fileIndex} = Inew;
        unAlignedImgs{fileIndex} = Iinit;
        e = Inew - U_hat{max_K} * status.kstatus{max_K}.w;
        outlierImgs{fileIndex} = e;

    if mod(fileIndex,SHOW_INTERVAL)==0 && ~OPTIONS.QUIET,
        % 1. initial unaligned image
        init_img = reshape(Iinit,imgSize(1),imgSize(2));
        axes(h_Initimg); imagesc(init_img);colormap gray;axis off;axis ij ;

        % 2. current aligned image
        alig_img = reshape( Inew ,imgSize(1),imgSize(2));
        axes(h_aligimg); imagesc(alig_img);colormap gray;axis off;axis ij ;

        % 3. lowrank reconstruction
        Uw_img = reshape(U_hat{max_K} * status.kstatus{max_K}.w, imgSize(1),imgSize(2));
        axes(h_training_uw); imagesc(Uw_img);colormap gray;axis off;axis ij ;

        % 4. outliers
      
        e_img = reshape(e,imgSize(1),imgSize(2)); 
        axes(h_outlier); imagesc(e_img);colormap gray;axis off;axis ij ;

    end
end
        
t_online = toc(t0);
fprintf('Online alignment %d frames: %.2f seconds\n', numImages,t_online);

end
