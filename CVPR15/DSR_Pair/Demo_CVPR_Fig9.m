close all;
dbstop if error

outdir = 'output';

figure; imshowpairflip(im, refim); savefig('noregistration.jpg', outdir); title('noregistration'); 
figure; imshowpairflip(newim3, refim); savefig('ssd.jpg', outdir);title('ssd');
figure; imshowpairflip(newim2, refim); savefig('RC.jpg', outdir);title('RC');
figure; imshowpairflip(newim1, refim); savefig('DTV.jpg', outdir); title('DTV');

% figure; imshowpairflip(warpAffine2(im_c,res3.X),refim,'blend'); title('ssd');


%  
% % [res2, newim2]=mirt2D_register(refim_dis,im, main, optim);
% I_RMSE3 = RMSE(newim3,refim);
% I_RMSE2 = RMSE(newim2,refim);
% % T_RMSE2 = RMSE(res2.X,T);
% % 
% I_RMSE1 = RMSE(newim1,refim);
% T_RMSE1 = RMSE(res1.X,T);
unregistered = imread('westconcordaerial.png');
figure, imshowflip(unregistered);  savefig('unregistered.jpg', outdir);
% text(size(unregistered,2),size(unregistered,1)+15, ...
%     'Image courtesy of mPower3/Emerge', ...
%     'FontSize',7,'HorizontalAlignment','right');
ortho = imread('westconcordorthophoto.png');
figure, imshowflip(ortho); savefig('westconcordorthophoto.jpg', outdir);
% text(size(ortho,2),size(ortho,1)+15, ...
%     'Image courtesy of Massachusetts Executive Office of Environmental Affairs', ...
%     'FontSize',7,'HorizontalAlignment','right');
load westconcordpoints;
input_points = movingPoints;
base_points = fixedPoints;
t_concord = cp2tform(input_points,base_points,'projective');
info = imfinfo('westconcordorthophoto.png');
registered = imtransform(unregistered,t_concord,...
                         'XData',[1 info.Width], 'YData',[1 info.Height]);
figure, imshowpairflip(registered,ortho,'blend'); savefig('MATLAB.jpg', outdir); savefig('MATres.jpg', outdir); title('MATLAB');



regsterRC = warpAffine2(im_c,res2);
regsterRC = regsterRC(1:size(ortho,1),1:size(ortho,2),:);
regsterDTV = warpAffine2(im_c,res1);
regsterDTV = regsterDTV(1:size(ortho,1),1:size(ortho,2),:);
figure; imshowpairflip(regsterRC,refim(1:size(ortho,1),1:size(ortho,2),:),'blend');  savefig('RCres.jpg', outdir);
figure; imshowpairflip(regsterDTV,refim(1:size(ortho,1),1:size(ortho,2),:),'blend'); savefig('DTVres.jpg', outdir);


y = 100:140; x = 240:280;
figure, imshowflip(unregistered(240:280,77:117,:)); savefig('ori_crop1.jpg', outdir);
figure, imshowpairflip(registered(x,y,:),ortho(x,y,:),'blend'); savefig('Mat_crop1.jpg', outdir); title('matlab');
figure; imshowpairflip(regsterRC(x,y,:),refim(x,y,:),'blend');  savefig('RC_crop1.jpg', outdir); title('RC');
figure; imshowpairflip(regsterDTV(x,y,:),refim(x,y,:),'blend'); savefig('DTV_crop1.jpg', outdir); title('DTV');

y = 280:320; x = 150:190;
figure, imshowflip(unregistered(120:160,258:298,:)); savefig('ori_crop2.jpg', outdir);
% figure, imshowflip(refim(x,y,:));
figure, imshowpairflip(registered(x,y,:),ortho(x,y,:),'blend'); savefig('Mat_crop2.jpg', outdir); title('matlab');
figure; imshowpairflip(regsterRC(x,y,:),refim(x,y,:),'blend'); savefig('RC_crop2.jpg', outdir); title('RC');
figure; imshowpairflip(regsterDTV(x,y,:),refim(x,y,:),'blend'); savefig('DTV_crop2.jpg', outdir); title('DTV');