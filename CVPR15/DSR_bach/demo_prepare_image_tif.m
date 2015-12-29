clear all; close all;

% fn = 'qb.jpg';
pre = 'Geo8';
% pre = 'iko4';
% pre = 'QB4';
% pre = 'SPOT1';
fn = [pre '.tif'];
im = imread(fn);

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
% subplot(n, m, 5); imshow(rgb2gray(im(:, :, 1:3)));
% subplot(n, m, 6); imshow(imadjust(rgb2gray(im), [], [], 1.3));
% subplot(n, m, 7); imshow(imadjust(rgb2gray(im), [], [], 0.7));
subplot(n, m, 8); imshow(rgb2gray(imlap(:, :, 1:3)));

if size(im, 3) == 3,
    imori = im;
    gray = rgb2gray(im(:,:,1:3));
    gray1 = rgb2gray(im(:,:,1:3));
    % gray2 = rgb2gray(im(:,:,2:4));
    % gray3 = rgb2gray(im(:,:,[1,3,4]));
    % gray4 = rgb2gray(im(:,:,[1,2,4]));
    imori(:,:,end+1) = gray1;
    % imori(:,:,end+1) = gray2;
    % imori(:,:,end+1) = gray3;
    % imori(:,:,end+1) = gray4;
    imori(:,:,end+1) = imadjust(gray, [], [], 1.2);
    imori(:,:,end+1) = imadjust(gray, [], [], 0.8);
    imori(:,:,end+1) = rgb2gray(imlap);
else
    imori = im;
    gray = rgb2gray(im(:,:,1:3));
    gray1 = rgb2gray(im(:,:,1:3));
    gray2 = rgb2gray(im(:,:,2:4));
    gray3 = rgb2gray(im(:,:,[1,3,4]));
    gray4 = rgb2gray(im(:,:,[1,2,4]));
    imori(:,:,end+1) = gray1;
    imori(:,:,end+1) = gray2;
%     imori(:,:,end+1) = gray3;
%     imori(:,:,end+1) = gray4;
end

save(['datasets/Satellite/' pre '.mat'], 'imori');

% figure; imhist(gray,64)

J = histeq(gray);
% figure; imhist(J)
% figure, imshow(J)

