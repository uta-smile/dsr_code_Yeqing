function res = showallimages(D,imagesize,cropsize)

[N,T] = size(D);
D = reshape(D,[imagesize(1) imagesize(2) 1 T]);

if nargin > 2
    D = D(cropsize+1:end-cropsize,cropsize+1:end-cropsize,1,:);
end


montage(D,'DisplayRange',[0 max(D(:))+eps]);