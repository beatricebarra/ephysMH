function [mu, sigma] = GMM_clustering(X, k)
    
    m = size(X, 1);
    % Randomly select k data points to serve as the means.
    indeces = randperm(m);
    mu = zeros(1, k);
    for (i = 1 : k)
        mu(i) = X(indeces(i));
    end
    
    % Use the overal variance ofhelp  the dataset as the initial variance for each cluster.
    sigma = ones(1, k) * sqrt(var(X));
    
    % Assign equal prior probabilities to each cluster.
    phi = ones(1, k) * (1 / k);
    
    %%===================================================
    % STEP 3: Run Expectation Maximization
    
    % Matrix to hold the probability that each data point belongs to each cluster.
    % One row per data point, one column per cluster.
    W = zeros(m, k);
    
    % Loop until convergence.
    for (iter = 1:10000)
        
        %fprintf('  EM Iteration %d\n', iter);
    
        %%===============================================
        % STEP 3a: Expectation
        %
        % Calculate the probability for each data point for each distribution.
        
        % Matrix to hold the pdf value for each every data point for every cluster.
        % One row per data point, one column per cluster.
        pdf = zeros(m, k);
        
        % For each cluster...
        for (j = 1 : k)
            
            % Evaluate the Gaussian for all data points for cluster 'j'.
            pdf(:, j) = gaussian1D(X, mu(j), sigma(j));
        end
        
        % Multiply each pdf value by the prior probability for each cluster.
        %    pdf  [m  x  k]
        %    phi  [1  x  k]   
        %  pdf_w  [m  x  k]
        pdf_w = bsxfun(@times, pdf, phi);
        
        % Divide the weighted probabilities by the sum of weighted probabilities for each cluster.
        %   sum(pdf_w, 2) -- sum over the clusters.
        W = bsxfun(@rdivide, pdf_w, sum(pdf_w, 2));
        
        %%===============================================
        % STEP 3b: Maximization
        %
        %Calculate the probability for each data point for each distribution.
    
        % Store the previous means so we can check for convergence.
        prevMu = mu;    
        
        % For each of the clusters...
        for (j = 1 : k)
        
            % Calculate the prior probability for cluster 'j'.
            phi(j) = mean(W(:, j));
            
            % Calculate the new mean for cluster 'j' by taking the weighted
            % average of *all* data points.
            mu(j) = weightedAverage(W(:, j), X);
        
            % Calculate the variance for cluster 'j' by taking the weighted
            % average of the squared differences from the mean for all data
            % points.
            variance = weightedAverage(W(:, j), (X - mu(j)).^2);
            
            % Calculate sigma by taking the square root of the variance.
            sigma(j) = sqrt(variance);
        end
        
        % Check for convergence.
        % Comparing floating point values for equality is generally a bad idea, but
        % it seems to be working fine.
        if (mu == prevMu)
            break
        end
    
    % End of Expectation Maximization loop.    
    end
end