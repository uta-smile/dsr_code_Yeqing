%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This is the demo code for


%%%% Contact: Chen Chen (cchen@mavs.uta.edu) and Junzhou Huang
%%%% (jzhuang@uta.edu), Department of Computer Science and Engineering,
%%%% University of Texas at Arlington



%%%% Partial of this code is from:
%%%% MIRT https://sites.google.com/site/myronenko/research/mirt
%%%% RASL http://perception.csl.illinois.edu/matrix-rank/rasl.html
%%%% SparseMRI http://www.eecs.berkeley.edu/~mlustig/Software.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, im_reg]=Register_DTV(refim, im, opts)
% refim: reference image
% im: source image
% opts: settings

% Output--
% T: the estimated transform
% im_reg: the image after registration


% checking for the other possible errors
if numel(size(im))~=numel(size(refim)), error('The dimensions of images are not the same.'); end;
if size(im,1)~=size(refim,1), error('The images must be of the same size.'); end;
if size(im,2)~=size(refim,2), error('The images must be of the same size.'); end;

% Original image size
imagesize=size(refim);
% dimen = length(imagesize);

% Create Pyramid
[fig1, fig2] = createPyramid(refim,im, opts.levels,opts.dimen);

if opts.dimen == 2
    Minitial = eye(3);
end
if opts.dimen == 3
    Minitial = eye(4);
end
if isfield(opts,'initialX')
    Minitial = opts.initialX;
end


% Go across all sub-levels
for level=1:opts.levels
    
    if opts.dimen == 2
        Minitial(1:2,3)=Minitial(1:2,3)*2;
    end
    if opts.dimen == 3
        Minitial(1:3,4)=Minitial(1:3,4)*2;
    end
    ima1 = fig1(opts.levels-level+1).im;
    ima2 = fig2(opts.levels-level+1).im;
    opts.T = Minitial;
    [M,result] = Registration(ima1,ima2,opts);
    
    Minitial = M;
end

% Prepair the output
T=M;

im_reg = result;
% because we have appended the initial images with the border of NaNs during the
% registration, now we want to remove that border and get the result of the
% initial image size
% im_reg=zeros(imagesize); [M,N,K]=size(result);
% im_reg(1:min(imagesize(1),M),1:min(imagesize(2),N),1:min(imagesize(3),K))=result(1:min(imagesize(1),M),1:min(imagesize(2),N),1:min(imagesize(3),K));

function [fig1, fig2] = createPyramid(im1,im2, level,dimen)
%CREATEPYRAMID Creates two n-level pyramid.
%   [fig1, fig2] = createPyramid(im1,im2,level) creates two pyramid of
%   depth defined by level. Returns two struct arrays. To access array
%   information use: fig1(1).im, fig1(2).im etc
%
%   See also IMRESIZE.
%   Copyright y@s
%   Date: Tuesday, Oct 22nd, 2002

% Assign lowest level of pyramid
fig1(1).im = im1;
fig2(1).im = im2;

% Loop to create pyramid
if dimen == 2
    for i=1: level-1
        fig1(1+i).im = imresize(fig1(i).im, [size(fig1(i).im,1)/2 size(fig1(i).im,2)/2], 'bilinear');
        fig2(1+i).im = imresize(fig2(i).im, [size(fig2(i).im,1)/2 size(fig2(i).im,2)/2], 'bilinear');
    end
end

if dimen == 3
    for i=1: level-1
        fig1(1+i).im = imresize3d(fig1(i).im, [],[size(fig1(i).im,1)/2 size(fig1(i).im,2)/2 size(fig1(i).im,3)/2], 'cubic');
        fig2(1+i).im = imresize3d(fig2(i).im, [],[size(fig2(i).im,1)/2 size(fig2(i).im,2)/2 size(fig2(i).im,3)/2], 'cubic');
    end
end
