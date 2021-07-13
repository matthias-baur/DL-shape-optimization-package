%% Gradient of eigenvalues wrt Fourier coefficients
% computes the gradient of eigenvalue lambda = k_crit^2 wrt Fourier coefficients

function d = eigval_gradient(a,b,param,k_crit, V)
    
    [r, ~, Gamma_col, Normal_col, p_source, ~, ~] = define_domain(a,b,param.t_col,param.p_beta);
    
    d = zeros(1,param.len_ab);
    lambda = k_crit^2;
    dn_eigf_temp = @(X,Y) dn_eigf(X,Y,k_crit, V,p_source,Normal_col); 
    
    % eigf_temp = @(X,Y) eigf(X,Y,k_crit, V,p_source);
    % dx = 10^(-4);
    % dn_eigf_temp = @(X,Y) dn_eigf_FD(X,Y,eigf_temp,Normal_col,dx);  % finite difference
    
    % evaluate formula for derivatives with quadrature
    % upper half with cos
    for j = 0:param.len_ab-1
        d(j+1) = param.dt_col*sum( (lambda - abs(dn_eigf_temp(Gamma_col(1,:),Gamma_col(2,:)).^2))  .*  (cos(j.*param.t_col)./r(param.t_col))  .*  sum(Gamma_col.*Normal_col,1) );
    end
    
    % lower half with sin
    for j = 0:param.len_ab-1
        d(j+1+param.len_ab) = param.dt_col*sum( (lambda - abs(dn_eigf_temp(Gamma_col(1,:),Gamma_col(2,:)).^2))  .*  (sin(j.*param.t_col)./r(param.t_col))  .*  sum(Gamma_col.*Normal_col,1) );
    end
    
end