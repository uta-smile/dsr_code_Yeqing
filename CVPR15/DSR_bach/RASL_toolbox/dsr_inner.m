function [R E deltaTau iter] = dsr_inner(D, J, alpha, beta, tol, maxIter,imgSize)

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
E = zeros( m, n);
R = D;
JdeltaTau = zeros( m, n);
obj_v = funv(R,E,D,JdeltaTau, alpha, beta,imgSize);
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

% d_norm = norm(D, 'fro');

iter = 0;
converged = false;
while ~converged
    iter = iter + 1;
    
     obj_v_old = obj_v;
    
     %update E  
     Z = D + JdeltaTau - R;
     E = wthresh(Z,'s',beta);
     %update R   
     Z = D + JdeltaTau - E;
     FZ = F(Z,imgSize);
     FZ(:,2:end) = wthresh(FZ(:,2:end),'s',alpha);
     R = FT(FZ,imgSize);
     
    %update deltaTau
    temp_T = R + E - D;
    
    for i = 1 : n
        deltaTau{i} = J{i}\temp_T(:,i);
        JdeltaTau(:, i) = J{i}*deltaTau{i} ;
    end
    
    obj_v = funv(R,E,D,JdeltaTau, alpha, beta,imgSize);

     
        
%     %update R,E
%     Z = D + JdeltaTau;
%     
%     W1 = F(R,imgSize);
%     W1 = (abs(W1)+eps).^(-1);
%        
%     W2 = abs(Z-R);
%     W2 = (abs(W2)+eps).^(-1);
%     
%     grad = FT(W1.*F(R,imgSize),imgSize) + lambda*W2.*R-lambda*W2.*Z;
%     
% 
%     for k=1:maxGradDescent
%         
% 
%         
%         R = R - gamma*grad;
%         E = D + JdeltaTau - R;
%         
%         obj_v = funv(R,E,lambda,imgSize);
%         
%         if obj_v>obj_v_old
%             R = R + gamma*grad;
%             gamma = gamma*beta;
% %             obj_v = obj_v_old;
%         else
%             obj_v_old = obj_v;
%         end
%         
%         
%     end
%     
%     %update deltaTau
%     temp_T = R + E - D;
%     
%     for i = 2 : n
%         deltaTau{i} = pinv(J{i})*temp_T(:,i);
%         JdeltaTau(:, i) = J{i}*deltaTau{i} ;
%     end
%     
%     obj_v = funv(R,E,lambda,imgSize);
    
    
    %     temp_T = D + dt_dual_matrix - E_dual + (1/mu)*Y;
    %     [U S V] = svd(temp_T, 'econ');
    %     diagS = diag(S);
    %     A_dual = U * diag(pos(diagS-1/mu)) * V';
    %
    %     temp_T = D + dt_dual_matrix - A_dual + (1/mu)*Y;
    %     E_dual = sign(temp_T) .* pos( abs(temp_T) - lambda/mu );
    %
    %     temp_T = D - E_dual - A_dual + (1/mu)*Y;
    %     for i = 1 : n
    %         dt_dual{i} =  - J{i}'*temp_T(:,i) ;
    %         dt_dual_matrix(:, i) = J{i}*dt_dual{i} ;
    %     end
    %
    %     Z = D + dt_dual_matrix - A_dual - E_dual;
    %     Y = Y + mu*Z;
    %
    %     obj_v = D(:)'*Y(:);
    %
    %     mu = mu*rho;
    %     stoppingCriterion = norm(Z, 'fro') / d_norm;
    stoppingCriterion = (obj_v_old - obj_v)/obj_v_old;
%     stoppingCriterion = 1;
    
    if mod( iter, DISPLAY_EVERY) == 0
        disp(['#Iteration ' num2str(iter)  ...
            ' objvalue ' num2str(obj_v) '  Stopping Criterion ' ...
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

