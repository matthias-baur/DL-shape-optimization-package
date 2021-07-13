%% Search for eigenvalues in interval
% searches for minima of f = log(|det(A(kappa))|) in interval

function [mins] = search_eigval_interval(f, interval, Nk, tol)

a=zeros(1,Nk+1);
kappa = linspace(interval(1), interval(2), Nk+1);

% evaluate function on interval in Nk equidistant points
for j = 0:Nk
    a(j+1) = f(kappa(j+1));
end

% plot f = log(|det(A)|) over kappa
% plot(kappa,a)
% hold on
% grid on

% find local maxima of function
k_local_max = [kappa(islocalmax(a)), interval(2)];

% search for local minima between subsequent local maxima and also last
% maximum and right interval boundary
mins = zeros(1, length(k_local_max)-1);
for j = 1:length(k_local_max)-1
    mins(j) = gr_min_search(f, k_local_max(j:j+1), tol);
end

% include leftmost interval if not beginning at zero
if interval(1) > 0 && ~isempty(k_local_max)
    mins = [gr_min_search(f, [interval(1) k_local_max(1)], tol), mins];
end

% exclude interval boundaries, not true local minima
mins = mins(mins-interval(1)>2*tol & interval(2)-mins>2*tol);

end