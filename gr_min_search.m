%% Golden ratio search for local minimum
% searches for a local minimizer of f:R->R in interval
% accuracy is set by tol

function [min] = gr_min_search(f, interval, tol)
    gr = (sqrt(5) + 1) / 2;
    a = interval(1);
    b = interval(2);
    
    c = b - (b - a) / gr;
    d = a + (b - a) / gr;
    while abs(b - a) > tol
        if f(c) < f(d)
            b = d;
        else
            a = c;
        end
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    end
    
    min = (b + a) / 2;
end
