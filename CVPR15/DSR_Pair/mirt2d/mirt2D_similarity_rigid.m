% mirt2D_similarity  The function computes the current similarity measure
% value and its dense gradients

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/

function [f,dx,imsmall]=mirt2D_similarity_rigid(im,refim,X,main)

% interpolate image and its gradients
% im_int=mirt2D_mexinterp(main.imsmall, Xx, Xy);
% clear Xx Xy tmp;
imsmall =  warpAffine2(im,X);

% imsmall=im_int(:,:,1); 
main.refimsmall = refim;
% Compute the similarity function value (f) and its gradient (dd)
switch lower(main.similarity)
   
   % sum of squared differences
   case 'ssd' 
        dd=imsmall-main.refimsmall;
        f=nansum(dd(:).^2)/2;
        
        ON = ~isnan(dd);
        ON = sum(ON(:));
        if ON == 0;
            f = inf;
        else
        f = f/ON;
        end
        
   % sum of absolute differences       
   case 'sad' 
        dd=imsmall-main.refimsmall;
        f=nansum(sqrt(dd(:).^2+1e-10));
        dd=dd./sqrt(dd.^2+1e-10); 
     
   % correlation coefficient     
   case 'cc' 
       
       %SJ=main.refimsmall-nansum(main.refimsmall(:))/numel(main.refimsmall);
       %SI=imsmall-nansum(imsmall(:))/numel(imsmall);
       mask=isnan(main.refimsmall+imsmall);
       main.refimsmall(mask)=nan;
       imsmall(mask)=nan;
       
       SJ=main.refimsmall-nanmean(main.refimsmall(:));
       SI=imsmall-nanmean(imsmall(:));
       
       
       a = nansum(imsmall(:).*SJ(:))/nansum(imsmall(:).*SI(:));
       f=-a*nansum(imsmall(:).*SJ(:));
       dd=-2*(a*SJ-a^2*SI);
       
    
   % Residual Complexity: A. Myronenko, X. Song: "Image Registration by
   % Minimization of Residual Complexity.", CVPR'09
   case 'rc' 
  
        rbig=imsmall-main.refimsmall;

        
        [y,x]=find_imagebox(rbig); r=rbig(y,x);
        r(isnan(r))=nanmean(r(:));

        Qr=mirt_dctn(r);
        Li=Qr.^2+main.alpha;

        f=0.5*sum(log(Li(:)/main.alpha));
        
        r=mirt_idctn(Qr./Li);
        dd=zeros(size(rbig)); 
        dd(y,x)=r;
        
                
        ON = ~isnan(rbig);
        ON = sum(ON(:));
        if ON == 0;
            f = inf;
        else
        f = f/ON;
        end
   % CD2 similarity measure: Cohen, B., Dinstein, I.: New maximum likelihood motion estimation schemes for
   % noisy ultrasound images. Pattern Recognition 35(2),2002
   case 'cd2' 
        
        f=(imsmall-main.refimsmall)/main.alpha;
        dd=2*tanh(f);
        f=2*nansum(log(cosh(f(:))));
        
   % MS similarity measure: Myronenko A., Song X., Sahn, D. J. "Maximum Likelihood Motion Estimation
   % in 3D Echocardiography through Non-rigid Registration in Spherical Coordinates.", FIMH 2009
   case 'ms' 
        f=(imsmall-main.refimsmall)/main.alpha;
        coshd2=cosh(f).^2;
        dd=tanh(f).*(2*coshd2+main.ro)./(coshd2-main.ro);
        f=nansum(1.5*log(coshd2(:)-main.ro)-0.5*log(coshd2(:)));
        
  % (minus) Mutual Information: Paul A. Viola "Alignment by Maximization of Mutual Information"   
  case 'mi' 
        % MI computation is somewhat more involved, so let's compute it in separate function       
        [f, dd]=mirt_MI(main.refimsmall,imsmall,main.MIbins);
        
  case 'dtv'
        %
        rbig=imsmall-main.refimsmall;
        [y,x]=find_imagebox(rbig); r=rbig(y,x);
        r(isnan(r))=nanmean(r(:));

        if size(r,1)<5
            f = inf;
            dx = zeros(size(rbig));
            return;
        end
        
        gradXY = main.TV*r;
        TV = (sum(abs(gradXY).^2,3) + 1e-15).^(1/2); 
        f = 0.5*sum(TV(:));

        dd=zeros(size(rbig)); 
        dd(y,x)=gTV(r,main);

        ON = ~isnan(rbig);
        ON = sum(ON(:));
               if ON == 0;
            f = inf;
        else
        f = f/ON;
        end
    otherwise
        error('Similarity measure is wrong. Supported values are RC, CC, MI, SAD, SSD, CD2, MS')
end

dx = dd;
% Multiply by interpolated image gradients
% ddx=dd.*im_int(:,:,2); ddx(isnan(ddx))=0;
% ddy=dd.*im_int(:,:,3); ddy(isnan(ddy))=0;

% This subfunctions finds the coordinates of the largest square
% within the image that has no NaNs (not affected by interpolation)
% It can be useful to ignore the border artifacts caused by interpolation,
% or when the actual image has some black border around it, that you don't
% want to take into account.
function [y,x]=find_imagebox(im)
[i,j]=find(~isnan(im)); 
n=4; % border size
y=min(i)+n:max(i)-n;
x=min(j)+n:max(j)-n;

function grad = gTV(x,main)
% compute gradient of TV operator

Dx = main.TV*x;
Dx2 = Dx.*conj(Dx);
G = (sum(Dx2,3) + 1e-15).^(1/2-1);
G = Dx.*repmat(G,[1,1,2]);
% G = Dx.*(Dx2 + 1e-15).^(1/2-1);
grad = main.TV'*G;

% function p=Grad1(u)
% % backward finite difference along dim 1
% %         u = reshape(u,n1,n2);
% p = [u(1,:);diff(u,1,1);];
% p = p(:);
% 
% 
% function q=Grad2(u)
% % backward finite difference along dim 2
% %         u = reshape(u,n1,n2);
% q = [u(:,1) diff(u,1,2)];
% q = q(:);
% 
% 
% function TVx=TV(x)
% % Total Variation norm of x, x is a n by n matrix
% grad_x = [Grad1(x) Grad2(x)];
% TV_eps = 0;
% if 1
%     pt_sqsum = sum(grad_x.*grad_x,2);
%     if TV_eps == 0; TVx = sum(sqrt(pt_sqsum)); else TVx = sum(sqrt(pt_sqsum+TV_eps)); end
% else
%     TVx = norm(grad_x(:),1);
% end
