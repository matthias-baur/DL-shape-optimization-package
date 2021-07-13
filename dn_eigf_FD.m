%% assumes given eigenfunction(x,y)

function z = dn_eigf_FD(x,y,eigf,normals,dx)

    z = zeros(size(x));
    for k = 1:length(x)
        z(k) = - ( eigf(x(k)-dx*normals(1,k),y(k)-dx*normals(2,k)) - eigf(x(k),y(k)) ) / dx;
    end

end