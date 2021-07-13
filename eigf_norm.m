%% L2 norm of eigenfunction
% approximates L2 norm of an approx. eigenfunction to kappa, defined by coefficients alpha
% and source points p_source, via simple quadrature

function L2norm = eigf_norm(kappa,alpha,p_source,dh)

    x=-1:dh:1;
    y=-1:dh:1;
    [X,Y] = meshgrid(x,y);
    Z = eigf(X,Y,kappa,alpha,p_source);

    L2norm = sqrt( dh^2*( norm(Z,'fro').^2 )); 
 
end