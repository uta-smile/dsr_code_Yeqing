function result = randpertube(im)

K = 40;
if nargin<3
    flag = 1;
end

[m,n,T]  = size(im);
result = im;

sigma = 3;



% distor = 0;
% result = im + 3*randn(size(im));
for k=1:K
    
    [Y,X] = meshgrid(1:n,1:m);
    mu_k = rand(2,1);
    mu_k(1) = round(mu_k(1)*m);
    mu_k(2) = round(mu_k(2)*n);
    
    Y = Y - mu_k(2);
    X = X - mu_k(1);
    
    z = X.^2 + Y.^2;
    z = -z/(2*sigma.^2);
    z = exp(z);
    
    z = z - mean(z(:));
    
    z = 3*z;
%     if flag == 1
%         distor = distor + z;
%     else
%         distor = distor - z;
%     end
    if mod(k,2) ==0
    result(:,:,1) = result(:,:,1) + z;
    else
        result(:,:,2) = result(:,:,2) + z;
%         result(:,:,1) = result(:,:,1) + z;
    end
end
if (flag ~= 1) & (K~=0)
distor = distor/K;
end

% result = result + 3*repmat(distor,[1,1,2]);
% result = result./(max(result(:)));
