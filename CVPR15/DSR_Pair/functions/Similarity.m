function [f,dx,imsmall]=Similarity(im,refim,X,opts)
% function f = ||T(s) - I||_TV
% Input:
%   im: source image I
%   refim: reference image R
%   X: transformation
%   main: settings


% Output:
%   f: function value
%   dx: the gradient of the function
%   imsmall:  % T(I)

imsmall =  warpAffine(im,X); % T(I)
% if opts.dimen ==2
%     tform = affine2d(X');
% end
% if opts.dimen ==3
%     tform = affine3d(X');
% end
% im = warpAffine2(im,X);
% imsmall = imwarp(im,tform);


opts.refimsmall = refim;
rbig=imsmall-opts.refimsmall;
% [y,x,z]=find_imagebox(rbig); r=rbig(y,x,z);
% r(isnan(r))=nanmean(r(:));
[nanpix]=find(isnan(rbig));
r = rbig;
r(nanpix) = zeros(size(nanpix));

if size(r,1)<5
    f = inf;
    dx = zeros(size(rbig));
    return;
end

if opts.dimen==3
    gradXY = TVOP3*r;
else
    gradXY = TVOP*r;
end
% TV = (abs(gradXY(:)).^2 + 1e-15).^(1/2);
% f = sum(TV); %L1 TV
TV = (sum(abs(gradXY).^2,opts.dimen+1) + 1e-15).^(1/2);
f = sum(TV(:));

% dd=zeros(size(rbig));
% dd(y,x)=gTV(r,opts);
dd = gTV(r,opts);

ON = ~isnan(rbig);
ON = sum(ON(:));
if ON == 0;
    f = inf;
else
    f = f/ON;
end


dx = dd;


% This subfunctions finds the coordinates of the largest square
% within the image that has no NaNs (not affected by interpolation)
% It can be useful to ignore the border artifacts caused by interpolation,
% or when the actual image has some black border around it, that you don't
% want to take into account.
function [y,x]=find_imagebox(im)
[i,j,k]=find(~isnan(im));
n=4; % border size
y=min(i)+n:max(i)-n;
x=min(j)+n:max(j)-n;
x=min(k)+n:max(k)-n;

function grad = gTV(x,opts)
% compute gradient of TV operator
if opts.dimen==3
Dx = TVOP3*x;
else
    Dx = TVOP*x;
end
Dx2 = Dx.*conj(Dx);
for k=1:size(Dx,opts.dimen+1)
    if opts.dimen==3
         G(:,:,:,k) = Dx(:,:,:,k).*(sum(Dx2,opts.dimen+1) + 1e-15).^(1/2-1);
    else
        G(:,:,k) = Dx(:,:,k).*(sum(Dx2,opts.dimen+1) + 1e-15).^(1/2-1);
    end

end
% G = Dx.*(Dx2 + 1e-15).^(1/2-1); %L1 TV
if opts.dimen==3
grad = TVOP3'*G;
else
   grad = TVOP'*G; 
    end
