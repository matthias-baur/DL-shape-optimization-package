%% Normal derivative of eigenfunction
% computes the normal derivative of an approx. eigenfunction to kappa, 
% defined by coefficients alpha and source points p_source, in points x,y on
% the boundary 

function z = dn_eigf(x,y,kappa,alpha,p_source,normals)

    % compute negative gradient
    gradient_u = zeros(2,length(x));
    dist_temp = sqrt(bsxfun(@minus, x, p_source(1,:)')'.^2+bsxfun(@minus, y, p_source(2,:)')'.^2);
    for k = 1:length(x)
        gradient_u(1,k) = (bsxfun(@minus, x(k), p_source(1,:)')' .* kappa .* (1i/4*besselh(1,kappa.* dist_temp(k,:)))./ dist_temp(k,:) ) * alpha ;
        gradient_u(2,k) = (bsxfun(@minus, y(k), p_source(2,:)')' .* kappa .* (1i/4*besselh(1,kappa.* dist_temp(k,:)))./ dist_temp(k,:) ) * alpha ;
    end
    
    % scalar product of gradients with normals
    z = - sum(gradient_u .* normals, 1);
end