function savemyfig(fn)

% disp(['skip ' fn])
% return

% saveas(gcf, fullfile('results', 'sate', fn));
F = getframe(gca);
Image = frame2im(F);
imwrite(Image, fullfile('results2', 'sate', fn))