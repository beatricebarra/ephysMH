function [rfModel, X_test, Y_test, Y_pred, nrmse] = emh_random_forest_regressor(X, Y, labels, label_symbols, scale)

% % Generate synthetic data (or replace with your dataset)
% rng(1); % Set random seed for reproducibility
% X = rand(100, 3) * 10;  % 100 samples, 3 features
% Y = 5 * X(:,1) - 3 * X(:,2) + 2 * X(:,3) + randn(100,1) * 2; % Noisy linear relationship

% Split into training and test sets
trainRatio = 0.8;
idx = randperm(length(Y));
trainIdx = idx(1:round(trainRatio * length(Y)));
testIdx = idx(round(trainRatio * length(Y)) + 1:end);

X_train = X(trainIdx, :);
Y_train = Y(trainIdx);
X_test = X(testIdx, :);
Y_test = Y(testIdx);
odor_test = labels(testIdx); 
% Train a Random Forest Regressor with 100 trees
numTrees = 1000;
rfModel = TreeBagger(numTrees, X_train, Y_train, 'Method', 'regression', 'OOBPrediction', 'on');

% Predict on test data
Y_pred = predict(rfModel, X_test);

% Evaluate performance
mse = mean((Y_test - Y_pred).^2);
rmse = sqrt(mean((Y_test - Y_pred).^2));
nrmse = rmse/(max([Y_test; Y_train]) - min([Y_test; Y_train]));
%mape = mean(abs((Y_test - Y_pred) ./ Y_test)) * 100;

fprintf('Mean Squared Error: %.3f\n', nrmse);

% Plot actual vs. predicted values

for i = 1:length(odor_test)
    iodor = odor_test(i); 
    mysymbol = label_symbols(iodor); 
    scatter(Y_test(i), Y_pred(i), 'filled', mysymbol, 'MarkerFaceColor', 'b');
end
hold on;

plot(Y_test, Y_test, 'r', 'LineWidth', 2); % Ideal prediction line
if strcmp(scale, 'log')
    %set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
end
xlabel('Actual Values');
ylabel('Predicted Values');
title('Random Forest Regression: Actual vs. Predicted');
grid on;

end