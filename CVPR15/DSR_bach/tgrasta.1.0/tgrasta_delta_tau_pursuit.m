%% t-GRASTA (TRANSFORMED GRASTA)
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


function [ w, e, dual, T_in_new, Inew, delta_norm ] = tgrasta_delta_tau_pursuit( U, I1, T_in, J, OPTS ,AlignPara, idx)
%% ADMM Solver for the Locally Linearized Problem


% transformation matrix to parameters
xi = projective_matrix_to_parameters(AlignPara.transformType,T_in) ;

% using QR to orthogonalize the Jacobian matrix
[Q, R] = qr(J,0) ;


[ w, e, delta_xi, dual, ~] = sparse_residual_pursuit(U(idx,:),  Q(idx,:),  I1(idx), OPTS);


xi = xi + R\delta_xi; %inv(R)*delta_xi ; %delta_xi;
T_in_new = parameters_to_projective_matrix(AlignPara.transformType, xi);% parameters to transformation matrix 

Inew = I1 + Q*delta_xi;

delta_norm = norm(delta_xi);

end

