function sol=RMSE(M,F)

idx1 = ~isnan(M);
idx2 = ~isnan(F);
idx = idx1 & idx2;


M=double(M(idx));
F=double(F(idx));
[n,m,d]=size(F);
D=(M(:,:,1:d)-F).^2;
sol=sqrt(sum(sum(sum(D)))/(n*m*d));
