


function J = image_Jaco3D(Ix, Iy, Iz, imgSize, transformType)

u   = vec(repmat(1:imgSize(2),imgSize(1),imgSize(3)));
v   = vec(repmat((1:imgSize(1))',imgSize(3),imgSize(2)));
z   = vec(repmat((1:imgSize(3)),imgSize(1)*imgSize(2),1)); 

if strcmp(transformType,'TRANSLATION'),
    J = [ Ix, Iy, Iz];
                    
elseif strcmp(transformType,'AFFINE'),
    J = [ Ix.*u,   Ix.*v,   Ix.*z, Ix,  Iy.*u,   Iy.*v, Iy.*z,  Iy, Iz.*u,   Iz.*v, Iz.*z, Iz];
                                
else
    error('Unrecognized transformation type in test_face_alignment.m');
end