function [y] = psychocurve(c, a, k, n)
    disp([a,k,n])
    for i = 1 : length(c)
        y(i) = a/(1 + exp(-k*(c(i)-n))); 
    end
    %y = a./(1 + exp(-k.*(c-n))); 
end