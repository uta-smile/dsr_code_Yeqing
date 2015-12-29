function [R E deltaTau iter] = dsr_inner_fast(D, J, alpha, beta, tol, maxIter,imgSize)

% Dec 2009
%
% by Chen Chen, 09, 2014 @UT Arlington cchen@mavs.uta.edus
% min  \alpha ||F*R||_1 + \\beta |E|_1  +  0.5*||R+E - D - J*deltaTau ||_2^2

[m n] = size(D);

if nargin < 7
    error('Too few arguments') ;
end

% if nargin < 4
%     tol = 1e-7;
% elseif tol == -1
%     tol = 1e-7;
% end
% 
% if nargin < 5
%     maxIter = 1000;
% elseif maxIter == -1
%     maxIter = 1000;
% end

DISPLAY_EVERY = 100 ;

% gamma = 1;
% beta = 0.001;
% alpha = 0.01;
% maxGradDescent = 200;

% initialize
% E = zeros( m, n);
E = zeros( m, n,'single');
R = D;
JdeltaTau = zeros( m, n,'single');
% obj_v = funv(R,E,D,JdeltaTau, alpha, beta,imgSize);
% Y = D;
% [U S V] = svd(D);
% Y = U * V';
% norm_two = norm(Y, 2);
% norm_inf = norm( Y(:), inf) / lambda;
% dual_norm = max(norm_two, norm_inf);
% Y = Y / dual_norm;
% norm_two = norm_two / dual_norm;
% norm_inf = norm_inf / dual_norm;
% obj_v = D(:)' * Y(:);

% A_dual = zeros( m, n);
% E_dual = zeros( m, n);
% dt_dual = cell(1,n) ;
% dt_dual_matrix = zeros(m, n) ;
% Jn = size(J{1}, 2) ; % the number of parameters of the transformation

% mu = 1.25/norm(D) ;
% rho = 1.25;

d_norm = norm(D, 'fro');
norm_Z = 1;
iter = 0;
converged = false;
while ~converged
    iter = iter + 1;
    norm_Z_old = norm_Z;
%      obj_v_old = obj_v;
    
     %update E  
     Z = D + JdeltaTau - R;
     E = wthresh(Z,'s',beta);
     %update R   
     Z = D + JdeltaTau - E;
     FZ = F(Z,imgSize);
%      FZ(:,2:end) = wthresh(FZ(:,2:end),'s',alpha);
     FZ = wthresh(FZ,'s',alpha);
     R = FT(FZ,imgSize);
     
    %update deltaTau
    temp_T = R + E - D;
    
    for i = 1 : n
        deltaTau{i} = J{i}'*temp_T(:,i);
        JdeltaTau(:, i) = J{i}*deltaTau{i} ;
    end
    
%     obj_v = funv(R,E,D,JdeltaTau, alpha, beta,imgSize);

    %     mu = mu*rho;
    %     stoppingCriterion = norm(Z, 'fro') / d_norm;
    Z = D + JdeltaTau - R - E;
    norm_Z = norm(Z);
    stoppingCriterion = abs((norm_Z - norm_Z_old)/norm_Z);
%     stoppingCriterion = norm(Z, 'fro') / d_norm;
%     stoppingCriterion = (obj_v_old - obj_v)/obj_v_old;
%     stoppingCriterion = 1;
    
    if mod( iter, DISPLAY_EVERY) == 0
        disp(['#Iteration ' num2str(iter)  ...
             '  Stopping Criterion ' ...
            num2str(stoppingCriterion)]);
    end
    
    
    if stoppingCriterion <= tol
        disp('RASL inner loop is converged at:');
        disp(['Iteration ' num2str(iter)  'Stopping Criterion'  ...
            num2str(stoppingCriterion)]) ;
        converged = true ;
    end
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;
    end
end


function out = funv(R,E,D,JdeltaTau, alpha, beta,imgSize)

v1 = F(R,imgSize);
v1 = alpha*sum(abs(v1(:)));
v2 = beta*sum(abs(E(:)));

v3 = norm(R+E-D-JdeltaTau,'fro');
v3 = 0.5*v3.^2;

out = v1+ v2+v3;

function out = F(X,imgSize)
[m n] = size(X);
X = reshape(X,[imgSize(1),imgSize(2),n]);
out = fft(X,[],3);
out = reshape(out,[m,n]);


function out = FT(X,imgSize)

[m n] = size(X);
X = reshape(X,[imgSize(1),imgSize(2),n]);
out = ifft(X,[],3);
out = reshape(out,[m,n]);

