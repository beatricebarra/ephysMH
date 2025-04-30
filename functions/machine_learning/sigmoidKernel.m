function K = sigmoidKernel(U, V)
    alpha = 0.01;  % scale
    c = 0;         % bias
    K = tanh(alpha * (U * V') + c);
end