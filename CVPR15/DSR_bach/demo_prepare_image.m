clear all; close all;

pre = 'geo10';
% pre = 'qb2';
fn = [ pre '.jpg'];
% fn = 'iko3.tif';
im = imread(fn);

stx = 201;
sty = 31;
len = 256;
im = im(stx:stx+len-1, sty:sty+len-1,:);

% h=fspecial('average',5);
h=fspecial('gaussian');
% h = fspecial('laplacian');
imlap = imfilter(im,h,'circular');

figure(1);
n = 2; m = 4;
subplot(n, m, 1); imshow(im(:, :, 1:3));
subplot(n, m, 2); imshow(im(:,:,1));
subplot(n, m, 3); imshow(im(:,:,2));
subplot(n, m, 4); imshow(im(:,:,3));
subplot(n, m, 5); imshow(rgb2gray(im(:, :, 1:3)));
subplot(n, m, 6); imshow(imadjust(rgb2gray(im), [], [], 1.2));
subplot(n, m, 7); imshow(imadjust(rgb2gray(im), [], [], 0.9));
subplot(n, m, 8); imshow(rgb2gray(imlap));

imori = im;
gray = rgb2gray(im);
imori(:,:,end+1) = mean(im, 3);
imori(:,:,end+1) = 0.8*im(:,:,1) + 0.1*im(:,:,1) + 0.1*im(:,:,3) ;
imori(:,:,end+1) = 0.1*im(:,:,1) + 0.8*im(:,:,1) + 0.1*im(:,:,3) ;
imori(:,:,end+1) = 0.1*im(:,:,1) + 0.1*im(:,:,1) + 0.8*im(:,:,3) ;
imori(:,:,end+1) = histeq(gray);
% imori(:,:,end+1) = imadjust(gray, [], [], 1.2);
% imori(:,:,end+1) = imadjust(gray, [], [], 0.9);
% imori(:,:,end+1) = rgb2gray(imlap);

save(['datasets/Satellite/' pre 'b.mat'], 'imori');
% imwrite(im, 'geo10b.jpg');

% figure; imhist(gray,64)

J = histeq(gray);
% figure; imhist(J)
% figure, imshow(J)

