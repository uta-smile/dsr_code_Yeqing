% t-GRASTA (TRANSFORMED GRASTA)
% Copyright (C) 2013-2014 by Jun He, Dejiao Zhang, and Laura Balzano
%  
%  This file is part of the t-GRASTA library.
%  It is provided without any warranty of fitness
%  for any purpose. You can redistribute this file
%  and/or modify it under the terms of the GNU
%  Lesser General Public License (LGPL) as published
%  by the Free Software Foundation, either version 3
%  of the License or (at your option) any later version.
%  (see http://www.opensource.org/licenses for more info)
%_________________________________________________________________________________________________________________________
% For t-GRASTA we align each frame with the canonical image size. We get
% the initial transformation--init_tau, by using similarity transformtion
% and the coordinates of two points choosen from the corresponding original 
% frame and canonical frame 
%
%         ______________________
%        |                      |                          
%        |                      |                   _______________
%        |    .           .     |   init tau       |(x1,y1)(x2,y2) |
%        |  (X1,Y1)    (X2,Y2)  |  ____________    |   .      .    |
%        |                      |                  |_______________|
%        |______________________|
%     
%     originalCoords = [X1 X2; ...          canonicalCoords = [ x1 x2; ...
%                         Y1 Y2 ];                                y1 y2  ];  
%
%  Note: For t-GRASTA we use multi_subspaces, the proper number of
%  subspaces is depended on the level of misalignment of certain dataset. 
%  The more severe misalignment the more subspaces we need.
%
%  We present three models of t-GRASTA, batch model, trained
%  online model and fully online model. You can test each of them. For both trained 
%  online model and fully online model, the superiority on efficiency and memory 
%  requirement become more and more obvious as the size of dataset become 
%  larger and larger.
%   
%  In this code, we contain five datasets. From "Gore", you can see t-GRASTA 
%  works well on video stabilization. For "Sidewalk","traffic", "Hall", 
%  t-GRASTA successfully accomplish the background/foreground separation even 
%  there exist different levels of video jitters. "Hall" contains severe 
%  misalignment, the jitters were manually simulated by us. Thus, it need 
%  more subspace to finish the separation well, approximately 20-25. "traffic",
%  "sidewalk" contain more misaligned frames which were caused by real video
%  jitters. As for them, we need less subspaces since the misalignment is 
%  not so severe as "Hall". 
%
%  "Digits" contains vary serious variation in each image and the dataset is 
%  small. Thus, it's difficult for us learn a pretty good subspace for 
%  either "fully online model" or "trained online model" by just subsampling a small
%  set of them. However, we can also efficiently align them by using the batch 
%  model of t-GRASTA. 
%  
% __________________________________________________________________________________________________________________________
%   
%   [1] Jun He, Dejiao Zhang, Laura Balzano, Tao Tao. Iterative Grassmannian Optimization for Robust Image Alignment.
%       preprint, http://arxiv.org/abs/1306.0404, June, 2013
%   [2] Jun He, Dejiao Zhang, Laura Balzano, Tao Tao. Iterative Online Subspace Learning for Robust Image Alignment, 
%       In IEEE Conference on Automatic Face and Gesture Recognition (FG), April 2013. 
%  
%  Author: Jun He, Dejiao Zhang, and Laura Balzano
%  Email:  hejun.zz@gmail.com, dejiaozhang@gmail.com, girasole@umich.edu
%  Project page: 
%           http://sites.google.com/site/hejunzz/tgrasta
%  Date:   July 04, 2013
% ___________________________________________________________________________________________________________________________
%
%  Disclaimer: 
%  Some basic image transform functions used in t-GRASTA depend on RASL, here we 
%  put those useful functions in the subdirectory 'RASL_toolbox_2010'.Interested users 
%  should refer to the authors' webpage for the latest version.
%       http://perception.csl.illinois.edu/matrix-rank/rasl.html
%   