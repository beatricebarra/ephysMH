function [net, X_test, Y_test, Y_pred, nrmse] = emh_perceptron_regressor(X, Y, labels, label_symbols, scale, NNeurons, trainRatio, transferFun)
    

    net = feedforwardnet(NNeurons);
    % Split dataset 
    % Split into training and test sets
    
    idx = randperm(length(Y));
    trainIdx = idx(1:round(trainRatio * length(Y)));
    testIdx = idx(round(trainRatio * length(Y)) + 1:end);
    
    X_train = X(trainIdx, :);
    Y_train = Y(trainIdx);
    X_test = X(testIdx, :);
    Y_test = Y(testIdx);
    odor_test = labels(testIdx); 
    % Change the output layer to a linear activation function for regression
    net.layers{end}.transferFcn = transferFun;
    
    % Train the network
    net = train(net, X_train', Y_train');
    
    % Test the network on new data
    
    Y_pred = net(X_test');
    % Evaluate performance
    mse = mean((Y_test - Y_pred).^2);
    rmse = sqrt(mean((Y_test - Y_pred).^2));
    nrmse = rmse/(max([Y_test; Y_train]) - min([Y_test; Y_train]));
    %mape = mean(abs((Y_test - Y_pred) ./ Y_test)) * 100;

    fprintf('Mean Squared Error: %.3f\n', nrmse);
    % Plor 
    for i = 1:length(odor_test)
        iodor = odor_test(i); 
        mysymbol = label_symbols(iodor); 
        scatter(Y_test(i), Y_pred(i), 'filled', mysymbol, 'MarkerFaceColor', 'b');
    end
    hold on;
    
    plot(Y_test, Y_test, 'r', 'LineWidth', 2); % Ideal prediction line
    if strcmp(scale, 'log')
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
    end
    axis equal
    xlabel('Actual Values');
    ylabel('Predicted Values');
    title('Perceptron prediction: Actual vs. Predicted');
    grid on;
end
