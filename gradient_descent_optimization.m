%% Gradient descent shape optimization
% Optimizes sums of eigenvalues by gradient descent method

% Domain initialization
a = A_pop(1,:); %zeros(1,len_ab);  %a_iter(17,:); %   
b = B_pop(1,:); %zeros(1,len_ab);  %b_iter(17,:); %   

%a(1) = 1/sqrt(pi);
%a(3) = 0.08;

% Initialize stepwise save variables
a_iter = a; b_iter = b;
d_iter = zeros(1,2*len_ab);
mult_iter = zeros(max_steps,n+1);
sum_n_iter = zeros(1,max_steps);

%% Gradient descent loop

while n_steps < max_steps
    tic
    % solve direct problem and compute gradient
    [sum_n, V, mult] = direct_problem(a,b,param,sum_length);
    
    sum_n_iter(n_steps+1) = sum_n; 
    mult_iter(n_steps+1,:) = mult;
    
    %if mult(1) > 1
    %   mult_tol = 0.7*mult_tol;
    %end
    
    d = eigval_gradient(a,b,param,mult(end),V(:,1));
    
    % add gradients of lower eigenvalues down to lowest in objective
    % function or lowest in [lambda_n-mult_tol , lambda_n]
    mult_number = max(mult(1),sum_length);
    for m = 2:mult_number
        d_temp = eigval_gradient(a,b,param,mult(end-(m-1)),V(:,m));
        d = d + d_temp;
    end
    d = d./mult_number;
    sum_n_iter(1:n_steps+1)
    
    % check on convergence
    if n_steps > 2 && abs(sum_n_iter(n_steps) - sum_n_iter(n_steps+1)) < eps
        break;
    end

    % optimize on line
    beta_max = 0.1*norm(a)/norm(d);
    beta_opt = gr_min_line_search(a,b,d,param, mult_number, [0 beta_max]);
    betas = [beta_opt beta_max]
    
    a = a - beta_opt*d(1:len_ab);
    b = b - beta_opt*d(len_ab+1:end);
    
    % rescale new coefficients
    [r, dr_dt] = define_r(a,b);
    [a,b] = rescale_domain(a,b,r,t_col,dt_col);
    
    % save new coefficients
    a_iter = [a_iter;a];
    b_iter = [b_iter;b];
    d_iter = [d_iter;d];
    
    n_steps = n_steps+1;
    mult_temp = mult;
    toc
    
end


