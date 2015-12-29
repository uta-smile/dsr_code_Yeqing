%% Aux function 
%
%%
function tgrasta_plot_bigpic( images, canonicalImageSize , plotTitle,bCell, bABS)

N = length(images);
xI = ceil(sqrt(N));
layout.xI = xI ;
layout.yI = xI ;
layout.gap = 0 ;
layout.gap2 = 0 ;

DIM = canonicalImageSize(1)*canonicalImageSize(2);
numImage = length(images);

if nargin < 5, bABS = 0; end

if bCell,
    imageMat = cell2mat(images);
    imageMat = reshape(imageMat,[DIM, numImage]);
else
    imageMat = images;
end

%% display

% layout
x = layout.xI ;
y = layout.yI ;
gap = layout.gap ;
gap2 = layout.gap2 ; % gap2 = gap/2;

container = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap); 
% white edges
bigpic = cell(x,y); % (xI*canonicalImageSize(1),yI*canonicalImageSize(2));

% D_unaligned
for i = 1:x
    for j = 1:y
        if y*(i-1)+j > numImage
            bigpic{i,j} = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap);
        else
            if bABS,
                a0 = abs(imageMat(:,y*(i-1)+j)) ; 
            else
                a0 = imageMat(:,y*(i-1)+j);
            end
            container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(a0, canonicalImageSize);
            bigpic{i,j} = container;
        end
    end
end
figure
imshow(cell2mat(bigpic),[],'DisplayRange',[0 max(max(imageMat))],'Border','loose')
title(plotTitle) ;
