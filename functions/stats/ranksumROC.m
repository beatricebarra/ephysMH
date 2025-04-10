function [X, Y, AUC, p] = ranksumROC(labels,scores)
    [X,Y,T,AUC] = perfcurve(labels,scores, 1); 
    
    
    n1 = sum(labels == 1); % Number of positive cases
    n0 = sum(labels == 0); % Number of negative cases
    
    % Compute Q1 and Q2
    Q1 = AUC / (2 - AUC);
    Q2 = (2 * AUC^2) / (1 + AUC);
    
    % Compute standard error (SE)
    SE_AUC = sqrt((AUC * (1 - AUC) + (n1 - 1) * (Q1 - AUC^2) + (n0 - 1) * (Q2 - AUC^2)) / (n0 * n1));
    
    % Compute Z-score
    Z = (AUC - 0.5) / SE_AUC;
    
    % Compute two-tailed p-value
    p = 2 * (1 - normcdf(abs(Z))); 
end