function result = add_distor(im,K,flag)

if nargin<3
    flag = 1;
end

[m,n,T]  = size(im);
result = im;

sigma = 30;



distor = 0;

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
    
    if flag == 1
        distor = distor + z;
    else
        distor = distor - z;
    end
    

end
if (flag ~= 1) & (K~=0)
distor = distor/K;
end

result = result + distor;
result = result./max(max(result));



