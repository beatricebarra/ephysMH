function my_linear_regression(x, y)
    % Put data in a table for fitlm
    T = table(x,y,'VariableNames',{'x','y'});
    
    % Fit linear model y ~ 1 + x
    mdl = fitlm(T,'y ~ x');
    
    % Display basic results
    disp(mdl)
    coef = mdl.Coefficients.Estimate;        % [Intercept; Slope]
    rsq = mdl.Rsquared.Ordinary;             % R-squared
    rmse = mdl.RMSE;                         % Root mean squared error
    
    fprintf('Intercept = %.4g, Slope = %.4g\n', coef(1), coef(2));
    fprintf('R^2 = %.4g, RMSE = %.4g\n', rsq, rmse);
    
    % Predictions and residuals
    yhat = predict(mdl, T);                  % fitted values
    res = y - yhat;                          % residuals
    
    % Optional: simple train/test split evaluation (70/30)
    cvFrac = 0.7;
    n = numel(x);
    idx = randperm(n);
    ntrain = round(cvFrac*n);
    trainIdx = idx(1:ntrain);
    testIdx  = idx(ntrain+1:end);
    
    mdl_cv = fitlm(x(trainIdx), y(trainIdx));    % fit on training set
    ytest_pred = predict(mdl_cv, x(testIdx));
    rmse_test = sqrt(mean((y(testIdx)-ytest_pred).^2));
    
    fprintf('Test RMSE (%.0f%% train): %.4g\n', 100*cvFrac, rmse_test);
    
    % Plot data, fit line, and residuals
    figure;
    subplot(2,1,1);
    scatter(x,y,15,'filled'); hold on;
    plot(x, predict(mdl, table(x,'VariableNames',{'x'})), 'r-', 'LineWidth',1.5);
    xlabel('x'); ylabel('y'); title('Data and Linear Fit');
    legend('Data','Fit','Location','best'); grid on;
    
    subplot(2,1,2);
    plot(x,res,'.');
    xlabel('x'); ylabel('Residual'); title('Residuals');
    yline(0,'r--'); grid on;
end