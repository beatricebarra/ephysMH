function my_random_forest(x, y)

% x: predictors (n×p) numeric matrix or column vector
% y: response (n×1) numeric vector
% Example:
% x = (0:0.1:10)'; y = 2.5*x + 1 + 0.5*randn(size(x));

% Ensure column orientation and sizes
x = double(x);
y = double(y(:));
assert(size(x,1)==numel(y), 'x and y must have the same number of rows.');

% If x is a single column vector, convert to table for clearer predictor names
if isvector(x) && size(x,2)==1
    T = table(x, y, 'VariableNames', {'x', 'y'});
    Xtrain = T(:, 'x');
else
    % If matrix, convert to table with generic names
    if ~istable(x)
        p = size(x,2);
        varNames = arrayfun(@(k) sprintf('x%d',k), 1:p, 'uni', false);
        Xtrain = array2table(x, 'VariableNames', varNames);
    else
        Xtrain = x; % already a table
    end
    Xtrain.y = y; % append response to keep consistent variable container if needed
    Xtrain = removevars(Xtrain, 'y'); % keep predictors only
end

% Train random forest (bagged ensemble of regression trees)
rng(0); % for reproducibility
numTrees = 200;                  % typical default
Mdl = fitrensemble(Xtrain, y, ...
    'Method','Bag', ...
    'NumLearningCycles',numTrees, ...
    'Learners','Tree'); 

% Predict on training data
yhat = predict(Mdl, Xtrain);

% Compute training RMSE
rmse_train = sqrt(mean((y - yhat).^2));

% k-fold cross-validated RMSE (k=5)
cvMdl = crossval(Mdl, 'KFold', 5);
yhat_cv = kfoldPredict(cvMdl);
rmse_cv = sqrt(mean((y - yhat_cv).^2));

% Display
fprintf('Train RMSE: %.4g\n', rmse_train);
fprintf('5-fold CV RMSE: %.4g\n', rmse_cv);

% Optional: predict for new data newX (matrix or table)
% ynew_pred = predict(Mdl, newX);
end