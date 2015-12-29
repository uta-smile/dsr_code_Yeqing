function [fx,fy,fz,ft] = computeDerivatives3(im1,im2)
% 
% function [fx,fy,fz] = computeDerivatives2(im1,im2)
%
% im1 and im2 are images
%
% [fx,fy,ft] are images, derivatives of the input images

filter = [0.03504 0.24878 0.43234 0.24878 0.03504];
dfilter = [0.10689 0.28461 0.0  -0.28461  -0.10689];

dx1 = conv3sep(im1,dfilter,filter,filter,'same');
dy1 = conv3sep(im1,filter,dfilter,filter,'same');
dz1 = conv3sep(im1,filter,filter,dfilter,'same');
blur1 = conv3sep(im1,filter,filter,filter,'same');
dx2 = conv3sep(im2,dfilter,filter,filter,'same');
dy2 = conv3sep(im2,filter,dfilter,filter,'same');
dz2 = conv3sep(im2,filter,filter,dfilter,'same');
blur2 = conv3sep(im2,filter,filter,filter,'same');

fx=(dx1+dx2)/2;
fy=(dy1+dy2)/2;
fz=(dz1+dz2)/2;
ft=(blur2-blur1);
return

function result = conv3sep(im,rowfilt,colfilt,depfilt,shape)
% CONV2SEP: Separable convolution using conv2.
% 
%      result=conv2sep(im,rowfilt,colfilt,shape)
%
%      im - input image.
%      rowfilt - 1d filter applied to the rows
%      colfilt - 1d filter applied to the cols
%      shape - 'full', 'same', or 'valid' (see doc for conv2).
%
% Example: foo=conv2sep(im,[1 4 6 4 1],[-1 0 1],'valid');


if ~exist('shape')
  shape='full';
end

rowfilt = reshape(rowfilt,[1,5,1]);
colfilt = reshape(colfilt,[5,1,1]);
depfilt = reshape(depfilt,[1,1,5]);

tmp = convn(im,rowfilt,shape);
result = convn(tmp,colfilt,shape);
result = convn(result,depfilt,shape);