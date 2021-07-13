%% Rescale domain
% rescales Fourier coefficients such that domain has area = 1

function [a1,b1] = rescale_domain(a,b,r,t_col,dt_col)

A_new = dt_col/2*(sum(r(t_col).^2));
a1 = a./sqrt(A_new);
b1 = b./sqrt(A_new);

end