% Aux function

function status = tgrasta_init_status(OPTIONS,DIM)
    status.curr_iter    = 0;
    status.last_mu      = 1;            
    status.step_scale   = 0;
    status.level        = 0;
    status.last_W       = zeros(OPTIONS.RANK, 1);
    status.last_gamma   = zeros(DIM,1);
    status.grad_ip      = 0;    
    status.OPTS.MAX_ITER = OPTIONS.ITER_MIN;
    status.OPTS.RHO      = OPTIONS.RHO;
    status.OPTS.TOL      = OPTIONS.INNER_TOL;
end
