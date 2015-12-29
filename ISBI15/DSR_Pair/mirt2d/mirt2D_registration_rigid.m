% MIRT2D_REGISTERATION  non-rigid registration
% of a single pair of 2D images at a given hierarchical level

% Copyright (C) 2007-2010 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/MIRT/
%
% This file is part of the Medical Image Registration Toolbox (MIRT).

function [X, im]=mirt2D_registration_rigid(refim,im, main, optim)

%% normalize the initial optimization step size
% compute the gradient of the similarity measure
% [Xx,Xy]=mirt2D_nodes2grid(X, main.F, main.okno);
% [f, ddx, ddy]=mirt2D_similarity(main, Xx, Xy);


% divide the initial optimization step size by the std of the gradient
% this somewhat normalizes the initial step size for different possible
% similarity measures used
% optim.gamma=optim.gamma/std([ddx(:); ddy(:)],0);
% clear ddx ddy Xx Xy

%% Start the main registration
% compute the objective function and its gradient
% [fe, T, im]=mirt2D_grad(X,  main);              % compute the similarity measure and its gradient
% [Xp, fr]=mirt2D_regsolve(X,T,main, optim, 1);   % compute the regularization term and update the transformation
% f=fe+fr;
% compute the value of the total objective function (similarity + regularization)
% main.refim = refim;
srim = im;

X = main.X;
% im = warpAffine2(im,X);

% delta_X = X;
% Compute the derivatives of the image
[Iu,Iv,ft]=computeDerivatives2(im,refim);
J =  image_Jaco(Iu(:), Iv(:), size(im), main.type, X);

% [f,ddx,ddy,imsmall]=mirt2D_similarity_rigid(im,main);
[f,dx,im]=mirt2D_similarity_rigid(im,refim,X,main);
%
% fchange=optim.fundif+1; % some big number
iter=0;
fchange  = inf;

if main.single
figure(22); imshowpair(im, refim);
pause(0.5);
end
% do while the relative function difference is below the threshold and
% the meximum number of iterations has not been reached
while (iter<optim.maxsteps) &&(abs(fchange)>optim.fundif)
    
    idx = ~isnan(J(:,1));   
    delta_X = eye(3);   
%     dx = dx(3:end-2,3:end-2);   
    idx2 = ~isnan(dx);
    idx =idx & idx2(:) ;
       
    A = J(idx,:);
    b = dx(idx);   
    %     [Q,R] = qr(A);
    grad = pinv(A)*b;
%     grad = A'*b;
    
    %     b = ft(idx);
    %     grad = A\b;
    converge = 0;iit = 0;optim.gamma = 1;
    while ~converge
        iit = iit + 1;
        
        if strcmp( main.type,'TRANSLATION'),
            delta_X(1:2,3) = delta_X(1:2,3) - optim.gamma*grad;
        elseif strcmp( main.type,'AFFINE'),
            %             grad = reshape(grad,[2,3]);
            delta_X(1,1:3) = delta_X(1,1:3) - optim.gamma*grad(1:3)';
            delta_X(2,1:3) = delta_X(2,1:3) - optim.gamma*grad(4:6)';
        end
        %         delta_X(1:2,1:3) = delta_X(1:2,1:3) - optim.gamma*reshape(grad,[2,3]);
        
        [fp,dpx,imb]=mirt2D_similarity_rigid(im,refim,delta_X,main);
        if (fp > f)
            
            if strcmp( main.type,'TRANSLATION'),
                delta_X(1:2,3) = delta_X(1:2,3) + optim.gamma*grad;
            elseif strcmp( main.type,'AFFINE'),
                %             delta_X(1:2,1:3) = delta_X(1:2,1:3) + optim.gamma*grad;
                %                 grad = reshape(grad,[2,3]);
                delta_X(1,1:3) = delta_X(1,1:3) + optim.gamma*grad(1:3)';
                delta_X(2,1:3) = delta_X(2,1:3) + optim.gamma*grad(4:6)';
                %                 delta_X(1:2,3) = delta_X(1:2,3) + optim.gamma*grad(:,1);
            end
            optim.gamma=optim.gamma*optim.anneal;
        else
            fnew = fp;
            break;
        end
        
        if iit>optim.maxsteps
            converge = 1;
            fnew = inf;
        end
        
    end
    
    X=delta_X*X;
    
    im = warpAffine2(srim,X);
    
    fchange=(fnew-f)/f;
    f=fnew; %im=imb;
    dx = dpx;
    
    if main.single
    figure(22); imshowpair(im, refim);
    pause(0.1);
    end
    
    [Iu,Iv,ft]=computeDerivatives2(im,refim);
    J =  image_Jaco(Iu(:), Iv(:), size(im), main.type, X);
    
    iter=iter+1;
    disp([iter, fnew]);
    
    if converge
        break;
    end
end

