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

function [ Unew, STATUSnew] = tgrasta_update( U_hat, STATUS, w, e, dual, y, OPTIONS, idx)
%  tgrasta_update is one gradient step along the geodesic of Grassmannian
%
%   Inputs:
%       U_hat: current estimated subspace
%       y:  current transformed vector
%       (w,e,dual): parameters estimated by ADMM solver    
%       STATUS:  the current status of t-GRASTA
%           last_mu: \mu_{t-1} for adaptive step rule
%           step_scale: the estimated step_scale constant for adaptive
%       step-size rule
%           last_gamma and last_w: previous gradient = last_gamma*last_w'
%
%   Outputs:
%       Unew:    the updated subspace
%       STATUSnew: the updated running status
%
%

QUIET               = 1;
LEVEL_FACTOR        = 2;
MAX_LEVEL           = 10;

MIN_MU              = 1;
MAX_MU              = 30; %15;
DEFAULT_MU_HIGH     = (MAX_MU-1)/2;
DEFAULT_MU_LOW      = MIN_MU + 2;

OPTS = STATUS.OPTS;

if isfield(OPTIONS,'ITER_MIN'),
    MIN_ITER = OPTIONS.ITER_MIN;
else
    MIN_ITER = 5;
end

if isfield(OPTIONS,'ITER_MAX'),
    ITER_MAX = OPTIONS.ITER_MAX;
else
    ITER_MAX = 20;
end

%update low rank matrix U for each scale
DIM = size(U_hat,1);

gamma_1 = dual + OPTIONS.RHO * (U_hat(idx,:) * w + e - y(idx)); % NOTE here: gamma_1 = dual also works!!!
UtDual = U_hat(idx,:)' * gamma_1;
gamma_2 = U_hat * UtDual;
gamma = zeros(DIM,1);
gamma(idx) = gamma_1;
gamma = gamma - gamma_2;

gamma_norm = norm(gamma);
W_norm     = norm(w);
sG         = gamma_norm * W_norm;

% Here we use the adaptive step-size rule for SGD. The adaptive can work
% well for both dynamic and static subspace tracking tasks

% 1. determine the step scale from the first observation
if ~STATUS.step_scale
    %     STATUS.STEP_SCALE_INIT   = 0.1*pi*(1+MIN_MU)/sG;
    STATUS.step_scale = 0.5*pi/sG;
    
    STATUS.constant_scale = 1*1e-3*pi/sG;
    
    if ~QUIET,
        fprintf('Constant C Estimation is: %.2e\n',STATUS.step_scale);
    end
end

% 2. inner product of previous grad and current grad
grad_ip = trace(STATUS.last_W * (STATUS.last_gamma' * gamma) * w');

%%% avoid inner product too large
normalization = norm(STATUS.last_gamma * STATUS.last_W','fro') * norm(gamma * w','fro');
if normalization == 0,
    grad_ip_normalization = 0;
else
    grad_ip_normalization = grad_ip/normalization;
end
%%%


% 3. if the two consecutive grad in the same direction, we take a larger
% step along the gradient direction, otherwise take a small step along the
% gradient direction
STATUS.last_mu = max(STATUS.last_mu + sigmoid(-grad_ip_normalization) , MIN_MU); %grad_ip_normalization

if OPTIONS.CONSTANT_STEP > 0,
    t = OPTIONS.CONSTANT_STEP; % STATUS.constant_scale * sG ; %
else
    % should not take a step larger than pi/2
    t = STATUS.step_scale * LEVEL_FACTOR^(-STATUS.level) * sG / (1+STATUS.last_mu);
    if t>=pi/3,
        t = pi/3;
    end
end
% fprintf('subspace update %.2e\n',t);


% Adjust the level(ScaleIndex)
bShrUpd = 0;
if STATUS.last_mu <= MIN_MU,
    if STATUS.level > 1,
        bShrUpd = 1;
        STATUS.level = STATUS.level - 1;
        
        if ~QUIET,
            fprintf('multi-level adaption - decreasing, t:%.2e, vectors: %d, level: %d\n',...
                t,STATUS.curr_iter, STATUS.level);
        end
        STATUS.curr_iter  = 0;
    end
    
    STATUS.last_mu   = DEFAULT_MU_LOW;
elseif STATUS.last_mu > MAX_MU
    if STATUS.level < MAX_LEVEL
        bShrUpd = 1;
        STATUS.level = STATUS.level + 1;
        
        if ~QUIET
            fprintf('multi-level adaption - increasing, t:%.2e, vectors: %d, level: %d\n',...
                t, STATUS.curr_iter,  STATUS.level);
        end
        STATUS.curr_iter  = 0;
        STATUS.last_mu    = DEFAULT_MU_HIGH;
    else
        STATUS.last_mu    = MAX_MU;
    end
end

if bShrUpd,
    if STATUS.level>=0 && STATUS.level <4,      % [0,2)
        OPTS.MAX_ITER = MIN_ITER;
    elseif STATUS.level>=4 && STATUS.level <7,      % [2,4)
        OPTS.MAX_ITER = min(MIN_ITER*2, ITER_MAX);
    elseif STATUS.level>=7 && STATUS.level <9,  % [4,6)
        OPTS.MAX_ITER = min(MIN_ITER*4, ITER_MAX);
    elseif STATUS.level>=9 && STATUS.level <10,  % [6,8)
        OPTS.MAX_ITER = min(MIN_ITER*8, ITER_MAX);
        %     elseif STATUS.level(ScaleIndex)>=8 && STATUS.level(ScaleIndex) <10, % [8,10)
        %         OPTS.OPTIONS.MAX_ITER = min(MIN_ITER*8, ITER_MAX);
        %     elseif STATUS.level(ScaleIndex)>=10 && STATUS.level(ScaleIndex) <12,% [10,12)
        %         OPTS.OPTIONS.MAX_ITER = min(MIN_ITER*10, ITER_MAX);
    else
        OPTS.MAX_ITER = ITER_MAX;               % [12,...)
    end
    if ~QUIET,    fprintf('Will use %d ADMM iterations\n',OPTS.MAX_ITER); end
end

% 4. update the gradient for further step size update
STATUS.last_gamma  = gamma;
STATUS.last_W      = w;

STATUS.grad_ip = grad_ip_normalization; %grad_ip_normalization just for debugging


% Take the gradient step along Grassmannian geodesic.
alpha = w/W_norm;
beta  = gamma/gamma_norm;
step  = (cos(t)-1)*U_hat*(alpha*alpha')  - sin(t)*beta*alpha';

U_hat = U_hat + step;

%%
e_full = zeros(DIM,1); e_full(idx) = e;
STATUS.e   = e_full;
STATUS.w     = w;
STATUS.ldual = dual;

STATUS.grasta_t = t;
STATUS.curr_iter = STATUS.curr_iter +1;

STATUS.OPTS = OPTS;
%%%%%%%%%%%%%%%%%%%%%%

Unew = U_hat;
STATUSnew = STATUS;
end



%% Function of Sigmoid

function fval = sigmoid(x)
FMIN = -1; FMAX = 1;
omega = 0.1; % 0.1

fval = FMIN + (FMAX - FMIN)/(1 - (FMAX/FMIN)*exp(-x/omega));
end