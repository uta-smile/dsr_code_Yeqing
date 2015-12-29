function f = reg_similarity(ref,im,type)



switch lower(type)
    
    % sum of squared differences
    case 'ssd'
        dd=ref-im;
        f=nansum(dd(:).^2)/2;
        
        % sum of absolute differences
    case 'sad'
        dd=ref-im;
        f=nansum(sqrt(dd(:).^2+1e-10));
        dd=dd./sqrt(dd.^2+1e-10);
        
        % correlation coefficient
    case 'cc'
        
        %SJ=main.refimsmall-nansum(main.refimsmall(:))/numel(main.refimsmall);
        %SI=imsmall-nansum(imsmall(:))/numel(imsmall);
        mask=isnan(ref+im);
        ref(mask)=nan;
        im(mask)=nan;
        
        SJ=ref-nanmean(ref(:));
        SI=im-nanmean(im(:));
        
        
        a = nansum(im(:).*SJ(:))/nansum(im(:).*SI(:));
        f=-a*nansum(im(:).*SJ(:));
        dd=-2*(a*SJ-a^2*SI);
        
        
        % Residual Complexity: A. Myronenko, X. Song: "Image Registration by
        % Minimization of Residual Complexity.", CVPR'09
    case 'rc'
        alpha = 0.05;
        rbig=im-ref;
        
        [y,x]=find_imagebox(rbig); r=rbig(y,x);
        r(isnan(r))=nanmean(r(:));
        
        Qr=mirt_dctn(r);
        Li=Qr.^2+alpha;
        
        f=0.5*sum(log(Li(:)/alpha));
        
        r=mirt_idctn(Qr./Li);
        dd=zeros(size(rbig));
        dd(y,x)=r;
        
        
        % CD2 similarity measure: Cohen, B., Dinstein, I.: New maximum likelihood motion estimation schemes for
        % noisy ultrasound images. Pattern Recognition 35(2),2002
    case 'cd2'
        alpha = 0.05;
        f=(im-ref)/alpha;
        dd=2*tanh(f);
        f=2*nansum(log(cosh(f(:))));
        
        % MS similarity measure: Myronenko A., Song X., Sahn, D. J. "Maximum Likelihood Motion Estimation
        % in 3D Echocardiography through Non-rigid Registration in Spherical Coordinates.", FIMH 2009
    case 'ms'
        alpha = 0.05;
        ro = 0.9;
        f=(im-ref)/alpha;
        coshd2=cosh(f).^2;
        dd=tanh(f).*(2*coshd2+ro)./(coshd2-ro);
        f=nansum(1.5*log(coshd2(:)-ro)-0.5*log(coshd2(:)));
        
        % (minus) Mutual Information: Paul A. Viola "Alignment by Maximization of Mutual Information"
    case 'mi'
        % MI computation is somewhat more involved, so let's compute it in separate function
        [f, dd]=mirt_MI(ref,im,64);
        
        
    case 'dtv'
        %
        [f] = 0.5*TV(ref-im);
        
    case 'ngc'
        [f] = NGC(ref,im);
    case 'l2g'
        [f] = 0.5*TV2(ref-im);
end

end


function p=Grad1(u)
% backward finite difference along dim 1
%         u = reshape(u,n1,n2);
p = [u(1,:);diff(u,1,1);];
p = p(:);
end

function q=Grad2(u)
% backward finite difference along dim 2
%         u = reshape(u,n1,n2);
q = [u(:,1) diff(u,1,2)];
q = q(:);
end

function TVx=TV(x)
% Total Variation norm of x, x is a n by n matrix
grad_x = [Grad1(x) Grad2(x)];
TV_eps = 0;
if 1
    pt_sqsum = sum(grad_x.*grad_x,2);
    if TV_eps == 0; TVx = sum(sqrt(pt_sqsum)); else TVx = sum(sqrt(pt_sqsum+TV_eps)); end
else
    TVx = norm(grad_x(:),1);
end
end

function [y,x]=find_imagebox(im)
[i,j]=find(~isnan(im));
n=4; % border size
y=min(i)+n:max(i)-n;
x=min(j)+n:max(j)-n;
end

function result=TV2(x)
% Total Variation norm of x, x is a n by n matrix
grad_x = [Grad1(x) Grad2(x)];
TV_eps = 0;
if 1
    pt_sqsum = sum(grad_x.*grad_x,2);
    if TV_eps == 0; result = sum((pt_sqsum)); else result = sum((pt_sqsum+TV_eps)); end
else
    result = norm(grad_x(:),1);
end
end

function result=NGC(I1,I2)
% Total Variation norm of x, x is a n by n matrix
grad1 = Grad1(I1)+ 1i*Grad2(I1);
grad2 = Grad1(I2)+ 1i* Grad2(I2);

grad1 = reshape(grad1,size(I1));
grad2 = reshape(grad2,size(I2));
R1 =  sqrt(Grad1(I1).^2 + Grad2(I1).^2);
R2 =  sqrt(Grad1(I2).^2 + Grad2(I2).^2);
R1 = reshape(R1,size(I1));
R2 = reshape(R2,size(I2));

% nm = ifft2(fft2(R1).*conj(fft2(R2)));
% gc = ifft2(fft2(grad1).*conj(fft2(grad2)));
% nm = conv2(R1,R2,'full');
% gc = conv2(grad1,conj(grad2),'full');
nm = R1.*R2;
gc = grad1.*conj(grad2);
result = sum(gc(:))./sum(nm(:));

% result = max(max(real(gc./nm)));
end


