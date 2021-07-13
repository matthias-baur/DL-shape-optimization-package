%% Random optimization
% slightly randomize coefficients and keep them if value of objective function 
% gets better

% domain initialization (take over final coeff. from gradient descent)
a = a_iter(end,:); 
b = b_iter(end,:); 
sum_n = sum_n_iter(size(a_iter,1));

% max deviations allowed 
dev_max = 0.025;

%% randomization loop

while n_steps < max_steps
    tic
    % randomize coefficients
    a_temp = (1+dev_max-2*dev_max*rand(1,len_ab)).*a;
    b_temp = (1+dev_max-2*dev_max*rand(1,len_ab)).*b;
    
    [sum_n_temp, V, mult_temp] = direct_problem(a_temp,b_temp,param,sum_length);
    
    % if sum is better, keep coefficients
    if sum_n_temp < sum_n
        sum_n = sum_n_temp
        mult = mult_temp;
        a = a_temp;
        b = b_temp;
        [r, dr_dt] = define_r(a,b);
        [a,b] = rescale_domain(a,b,r,t_col,dt_col);
    end
    n_steps = n_steps+1;
    toc
end