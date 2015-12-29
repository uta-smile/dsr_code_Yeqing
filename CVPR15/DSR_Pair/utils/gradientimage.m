function out = gradientimage(im)

% [m,n] = size(im);

grad = TVOP*im;


out = sqrt(sum(grad.^2,3));
