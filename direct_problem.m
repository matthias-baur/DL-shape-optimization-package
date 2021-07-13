%% finds eigenfrequencies for given domain 
% last argument can be used to give only n-th eigenfreq

function [sum_n, V, mult] = direct_problem(a,b,param,sum_length)

% rescale domain to area=1
[r, ~] = define_r(a,b);
[ca,cb] = rescale_domain(a,b,r,param.t_col,param.dt_col);
[~, ~, Gamma_col, ~, p_source, ~, ~] = define_domain(ca,cb,param.t_col,param.p_beta);

% Build matrix function of kappa (imposing boundary conditions)
diff_matrix_x = bsxfun(@minus, Gamma_col(1,:), p_source(1,:)')';
diff_matrix_y = bsxfun(@minus, Gamma_col(2,:), p_source(2,:)')';
A = @(kappa) 0.5*sqrt(length(Gamma_col))*1i/4*besselh(0,kappa.* sqrt(diff_matrix_x.^2+diff_matrix_y.^2));

log_detA = @(kappa) log(abs(det(A(kappa))));

% find eigenfrequencies by searching log_detA for minima
interval = [0 param.interval_end];  % search interval for eig.frq
Nk = 400;                           % no. of 
tol = 10^(-8);                      % eig.frq accuracy tolerance

k_crit = search_eigval_interval(log_detA, interval, Nk, tol);

% search iteratively on finer mesh if necessary
dk = 0.1*(interval(2)-interval(1));
while length(k_crit)<param.n
        if dk < tol
            error = "could not find k_crit(n), dk < tol"
            sum_n = k_crit(param.n-1);
            break
        end
        % search around already found minima 
        for m = 2:length(k_crit)
            interval = [k_crit(m)-dk, min(k_crit(m)+dk,param.interval_end)];
            k_crit_temp = search_eigval_interval(log_detA, interval, Nk, tol);
            if length(k_crit_temp) > 1
                % filter out already found k_crit
                dist_to_k_crit = min(abs(k_crit_temp-k_crit'),[],1);
                k_crit_temp=k_crit_temp(dist_to_k_crit>tol & k_crit_temp > 4.2);
                % insert new k_crit
                k_crit = sort([k_crit, k_crit_temp]);
            end
        end
        dk=0.1*dk;               
end   

% if not enough eigenfrequencies found, fill up k_crit with param.interval_end
if length(k_crit)<param.n
    k_crit = [k_crit param.interval_end*ones(1,n-length(k_crit))];
end

% save multiplicity of n-th eigenvalue and all eig.frequencies in mult
mult = [sum(abs(k_crit.^2-k_crit(param.n).^2) < param.mult_tol), k_crit(1:param.n)];

% compute coefficients for L2-normalized eigenfunctions
dh = 0.01;
V = zeros(length(param.t_col),mult(1));    
[V_temp,~] = eigs(A(k_crit(param.n)),1,'smallestabs');
L2norm = eigf_norm(k_crit(param.n),V_temp,p_source,dh);
V(:,1) = V_temp./L2norm;
for m = 2:max(mult(1), sum_length)
    [V_temp,~] = eigs(A(mult(end-(m-1))),1,'smallestabs');
    L2norm = eigf_norm(mult(end-(m-1)),V_temp,p_source,dh);
    V(:,m) = V_temp./L2norm;
end

% compute sum of eigenvalues (value of objective function)
sum_n = sum(mult(param.n+1-(sum_length-1):param.n+1).^2); 


% Polya conjecture check. Not expected to reach interior of loop if
% conjecture is true. Can happen numerically, when bad parameters are
% chosen and the eigenfrequencies are very off from their true values.
if sum( k_crit(1:param.n) < sqrt(4*pi.*(1:param.n)) ) > 0
    disp("----------------------------------")
    disp("!!!Polya's conjecture violation!!!")
    disp("----------------------------------")
    a
    b
    [k_crit(1:param.n); sqrt(4*pi.*(1:param.n))]
end


end
