function imshowpairflip(im1, im2, arg)

if exist('arg', 'var'),
    imshowpair(flip(im1), flip(im2), arg)
else 
    imshowpair(flip(im1), flip(im2))
end