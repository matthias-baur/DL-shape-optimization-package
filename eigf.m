%% Eigenfunction
% evaluates an approx. eigenfunction to kappa, defined by coefficients alpha
% and source points p_source, in points X, Y (meshgrid)
function Z = eigf(X,Y,kappa,alpha,p_source)

    Z = zeros(size(X));
    for k = 1:size(X,1)
        Z(k,:) = (1i/4*besselh(0,kappa.* sqrt(bsxfun(@minus, X(k,:), p_source(1,:)')'.^2+bsxfun(@minus, Y(k,:), p_source(2,:)')'.^2))) * alpha;
    end

end