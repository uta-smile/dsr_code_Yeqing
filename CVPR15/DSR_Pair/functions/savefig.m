function savefig(name, od)

if ~exist('od', 'var'), 
    od = '.';
end
fn = fullfile(od, name);
saveas(gcf, fn);

% crop the interest area
im = imread(fn);
[m, n, k] = size(im);
[x1, y1] = deal(1);
[x2, y2] = deal(m,n);
for i=1:k,
    y1 = max(y1, find(im(m/2, :, i)<248, 1, 'first'));
    y2 = min(y2, find(im(m/2, :, i)<248, 1, 'last'));
    
    x1 = max(x1, find(im(:, n/2, i)<248, 1, 'first'));
    x2 = min(x2, find(im(:, n/2, i)<248, 1, 'last'));
end

im = im(x1:x2, y1:y2, :);
% figure; imshow(im)
imwrite(im, fn);

