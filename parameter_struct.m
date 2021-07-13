%% Parameter struct 
% creates a struct that contains the most important parameters
function param = parameter_struct(n,sum_length,interval_end,len_ab,N_col,dt_col,t_col,p_beta,beta_opt,beta_rel,mult_tol)
    param.n = n;
    param.sum_length = sum_length;
    param.interval_end = interval_end;
    param.len_ab = len_ab;
    param.N_col = N_col;
    param.dt_col = dt_col;
    param.t_col = t_col;
    param.p_beta = p_beta;
    param.beta_opt = beta_opt;
    param.beta_rel = beta_rel;
    param.mult_tol = mult_tol;
end