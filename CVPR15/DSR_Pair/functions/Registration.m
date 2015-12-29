function [X, im]=Registration(refim,im, opts)
% image registration with DTV
% min f = ||T(s) - I||_TV



%% Start the main registration
srim = im;

X = opts.T;
% if opts.dimen ==2
%     tform = affine2d(X');
% end
% if opts.dimen ==3
%     tform = affine3d(X');
% end
im = warpAffine(im,X);
% im = imwarp(im,tform);

% Compute the derivatives of the image
    if opts.dimen ==2
    [Iu,Iv,ft]=computeDerivatives2(im,refim); %%% need to be extended to 3D
    J =  image_Jaco(Iu(:), Iv(:), size(im), opts.type, X); %%% need to be extended to 3D
    end
    
    if opts.dimen ==3
        [Iu,Iv,fz,ft]=computeDerivatives3(im,refim); %%% need to be extended to 3D
        J =  image_Jaco3D(Iu(:), Iv(:), fz(:), size(im), opts.type); %%% need to be extended to 3D
    end

[f,dx,im]=Similarity(im,refim,X,opts);  %%% need to be extended to 3D

iter=0;
fchange  = inf;

if opts.display &&  opts.dimen ==2
    figure(22); imshowpair(im, refim);
    pause(0.5);
end
% do while the relative function difference is below the threshold and
% the meximum number of iterations has not been reached
while (iter<opts.maxsteps) &&(abs(fchange)>opts.fundif)
    
    idx = ~isnan(J(:,1));
    
    
    delta_X = eye(opts.dimen+1);
%     dx = dx(3:end-2,3:end-2);
    idx2 = ~isnan(dx);
    idx =idx & idx2(:) ;
    
    A = J(idx,:);
    b = dx(idx);
    grad = pinv(A)*b;
    
    converge = 0;iit = 0;opts.gamma = 1;
    while ~converge
        iit = iit + 1;
        if opts.dimen ==2
            if strcmp( opts.type,'TRANSLATION'),
                delta_X(1:2,3) = delta_X(1:2,3) - opts.gamma*grad;
            elseif strcmp( opts.type,'AFFINE'),
                delta_X(1,1:3) = delta_X(1,1:3) - opts.gamma*grad(1:3)';
                delta_X(2,1:3) = delta_X(2,1:3) - opts.gamma*grad(4:6)';
            end
        end
        if opts.dimen ==3
            if strcmp( opts.type,'TRANSLATION'),
                delta_X(1:3,4) = delta_X(1:3,4) - opts.gamma*grad;
            elseif strcmp( opts.type,'AFFINE'),
                delta_X(1,1:4) = delta_X(1,1:4) - opts.gamma*grad(1:4)';
                delta_X(2,1:4) = delta_X(2,1:4) - opts.gamma*grad(5:8)';
                delta_X(3,1:4) = delta_X(3,1:4) - opts.gamma*grad(9:12)';
            end
        end
        
        
        [fp,dpx,imb]=Similarity(im,refim,delta_X,opts);
        
        
        if (fp > f) % backtracking
            if opts.dimen ==2
                if strcmp( opts.type,'TRANSLATION'),
                    delta_X(1:2,3) = delta_X(1:2,3) + opts.gamma*grad;
                elseif strcmp(opts.type,'AFFINE'),
                    delta_X(1,1:3) = delta_X(1,1:3) + opts.gamma*grad(1:3)';
                    delta_X(2,1:3) = delta_X(2,1:3) + opts.gamma*grad(4:6)';
                end
            end
            
        if opts.dimen ==3
            if strcmp( opts.type,'TRANSLATION'),
                delta_X(1:3,4) = delta_X(1:3,4) + opts.gamma*grad;
            elseif strcmp( opts.type,'AFFINE'),
                delta_X(1,1:4) = delta_X(1,1:4) + opts.gamma*grad(1:4)';
                delta_X(2,1:4) = delta_X(2,1:4) + opts.gamma*grad(5:8)';
                delta_X(3,1:4) = delta_X(3,1:4) + opts.gamma*grad(9:12)';
            end
        end
            
            opts.gamma=opts.gamma*opts.anneal;
        else
            fnew = fp;
            break;
        end
        
        if iit>opts.maxsteps
            converge = 1;
            fnew = inf;
        end
        
    end
    
    X=delta_X*X;
    
%     if opts.dimen ==2
%         tform = affine2d(X');
%     end
%     if opts.dimen ==3
%         tform = affine3d(X');
%     end
%     im = imwarp(srim,tform);
        im = warpAffine(srim,X);
    
    fchange=(fnew-f)/f;
    f=fnew;
    dx = dpx;
    
    if opts.display &&  opts.dimen ==2
        figure(22); imshowpair(imb, refim);
        pause(0.1);
    end
    
    if opts.dimen ==2
    [Iu,Iv,ft]=computeDerivatives2(im,refim); %%% need to be extended to 3D
    J =  image_Jaco(Iu(:), Iv(:), size(im), opts.type, X); %%% need to be extended to 3D
    end
    
    if opts.dimen ==3
        [Iu,Iv,fz,ft]=computeDerivatives3(im,refim); %%% need to be extended to 3D
        J =  image_Jaco3D(Iu(:), Iv(:), fz(:), size(im), opts.type); %%% need to be extended to 3D
    end
    
    iter=iter+1;
    
    if converge
        break;
    end
end

