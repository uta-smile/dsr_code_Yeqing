%% Aux function for making the demo video

function make_video( video_name, video_matrix_fg,video_matrix_bg, aligned_matrix, video_matrix,vInfo)

frame_count = size(video_matrix_fg,1);
rows = vInfo.rows; cols = vInfo.cols;
vidObj = VideoWriter([video_name '.avi']);
open(vidObj);

figure(1); 

hvideo = subplot(2,2,1);set(gca,'nextplot','replacechildren');title('I');
halign = subplot(2,2,2);set(gca,'nextplot','replacechildren');title('I o tau');
hbg = subplot(2,2,3);set(gca,'nextplot','replacechildren');title('Uw');
hfg = subplot(2,2,4);set(gca,'nextplot','replacechildren');title('e');



colormap gray;axis off;
clims = [-0.8 0.8];
for i=1:frame_count,
    
    img = reshape(video_matrix_fg{i},rows,cols);
%     mmax = max(abs(img(:))); 
%     if  mmax < 0.2,
%         img(1,1) = 1;
%     end
    
    axes(hfg); imagesc(50*abs(img),clims); colormap gray;axis off; axis ij ; 
    
    img = reshape(video_matrix_bg{i},rows,cols);
    axes(hbg); imagesc((img)); colormap gray;axis off; axis ij ; 
    
    img = reshape(aligned_matrix{i},rows,cols);
    axes(halign); imagesc((img)); colormap gray;axis off; axis ij ; 
    
    img = reshape(video_matrix{i},rows,cols);
    axes(hvideo); imagesc((img)); colormap gray;axis off; axis ij ;  

    % Write each frame to the file.
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end

% Close the file
close(vidObj);
fprintf('Video have been written into %s video file\n',video_name);


end

