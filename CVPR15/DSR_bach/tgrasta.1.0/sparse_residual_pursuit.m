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
%  (see http://www.opensource.org/licenses for more info

function [ w, e, delt_tau, y, numInnerIter ] = sparse_residual_pursuit( A, J, v_hat, OPTS)
% Solve sparse residual pursuit (SRP) via ADMM {S.Boyd 2011}
%
% [ w, e, delt_tau, numIterInnerEach_1 ] = sparse_residual_pursuit( A, v_hat, J, OPTS)
% 
% Solves the following problem via ADMM:
% 
%   minimize     ||e||_1
%   subject to   v_hat + J*delt_tau - Aw -e = 0
%   y is the dual vector . 
%
%
%Batch image alignment: Alig_Relig inner loop
%
% input: A                       --- the matrix with vectors of image of only one subject  
%        J                       --- the Jacobian matrix
%        v_hat                   --- the normalize vector of testing image aligned by the previous tau      
%        OPTIONS                   --- structure for tuning the behavior of SRP
%
% output: w                      ---
%         e                      --- 
%         delt_tau               ---  
%         y                      --- dual value
%         numIterInner           --- total number of inner loop iterations
%
% RHO   : the augmented Lagrangian parameter, default 1.8
% MAX_ITER: max iteration for SRP problem
%   


%% Global constants and defaults

if isfield(OPTS,'RHO'),
    rho = OPTS.RHO;
else
    rho = 1.8;
end

if isfield(OPTS,'MAX_ITER'),
    MAX_ITER = OPTS.MAX_ITER;
else
    MAX_ITER = 100;
end


if isfield(OPTS,'TOL'),
    tol = OPTS.TOL;
else
    tol = 1e-6;
end

%% Data preprocessing

[m, n] = size(A);

% precompute static variables for w-update 
P = (A'*A) \ (A');

%precompute static variables for delt_tau update
F = (J'*J) \ (J');

mu = rho; %1.25/norm(v_hat);

v_norm = norm(v_hat);

%prepare parameters for termination checks
e = zeros(m,1);
y = zeros(m,1);
w = zeros(n,1);

iter = 0;
converged = false;

while ~converged && iter < MAX_ITER,
    iter = iter + 1;
    %delt_tau update
    delt_tau = F * (A*w + e - v_hat + y/mu);
    
    % w update
    w = P * (v_hat + J*delt_tau - e - y/mu);
    
    % e update
    Aw_hat = A*w;
    e = shrinkage( v_hat + J*delt_tau - Aw_hat - y/mu, 1/mu);
    
    % h update
    
    h = Aw_hat + e -v_hat - J*delt_tau;
    
    y = y + h*mu;
    
    mu = mu * 2; % for acclerating convergence
    
    stoppingCriterion = norm(h) / v_norm;
    
    if stoppingCriterion < tol,
        converged  = 1;
    end
           
end

numInnerIter= iter ;

end



function y = shrinkage(a, kappa)
    y = max(0, a-kappa) - max(0, -a-kappa);
end