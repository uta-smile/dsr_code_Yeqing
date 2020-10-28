clear all; close all;

DATA_PATH = '/Users/yeqing/Documents/Data/multipie/session01';
num_subjects = 3;

allboxes = {};

for sub = 1:num_subjects
    %% define images' path
    
    subject = num2str(sub,'%03d');
    
    path = fullfile(DATA_PATH, subject, '01', '05_1');
    Files = dir(fullfile(path,'*.png'));
    LengthFiles = length(Files);
    for i=1:LengthFiles
        filename = Files(i).name;
        filename = fullfile(path, filename);
        im = double(imread(filename));
        
        im = mean(im,3);
        im = imresize(im,[240 320]);
        I0(:,:,i) = im;
        
    end
    
    faceDetector = vision.CascadeObjectDetector;
    %bboxes = step(faceDetector, I0(:,:,1)/255); % Before Matlab 2016
    bboxes = faceDetector(I0(:,:,1)/255); % After MatLab 2016
    allboxes{end+1} = bboxes; 
    save allboxes allboxes
end

disp('Done')