%% Golden ratio search for local maximum
% searches for a local maximizer of f:R->R in interval
% accuracy is set by tol

function [max] = gr_max_search(f, interval, tol)
    gr = (sqrt(5) + 1) / 2;
    a = interval(1);
    b = interval(2);
    
    c = b - (b - a) / gr;
    d = a + (b - a) / gr;
    while abs(b - a) > tol
        if f(c) > f(d)
            b = d;
        else
            a = c;
        end
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    end
    
    max = (b + a) / 2;
end
