function [CVEV, MSE, Y_pred_all] = myRidge(X_train_all, Y_train_all, X_test_all, Y_test_all, lambda)
    % Ridge regression for multiple behaviors (loco, wsk, state)
    % Handles NaNs, independent X inputs, and outputs metrics for each behavior
    % 
    % Inputs:
    %   - X_train_all: training data for each behaviour {X_train_loco, X_train_wsk, X_train_state}
    %   - Y_train_all: training labels {Y_train_loco, Y_train_wsk, Y_train_state}
    %   - X_test_all: testing data for each behaviour {X_test_loco, X_test_wsk, X_test_state}
    %   - Y_test_all: testing labels {Y_test_loco, Y_test_wsk, Y_test_state}
    %   - lambda: Ridge regression regularisation parameter
    % 
    % Outputs:
    %   - CVEV: explained variances for each behavior
    %   - MSE: mean squared errors for each behavior
    %   - Y_pred_all: Cell array of predicted values for each behavior

    num_behaviors = length(Y_train_all);
    CVEV = zeros(1, num_behaviors);
    MSE = zeros(1, num_behaviors);
    Y_pred_all = cell(1, num_behaviors);

    for i = 1:num_behaviors

        X_train = X_train_all{i};
        Y_train = Y_train_all{i};
        X_test = X_test_all{i};
        Y_test = Y_test_all{i};

        % Remove NaNs
        valid_train = ~any(isnan(X_train), 2) & ~isnan(Y_train);
        valid_test = ~any(isnan(X_test), 2) & ~isnan(Y_test);
        X_train = X_train(valid_train, :);
        Y_train = Y_train(valid_train);
        X_test = X_test(valid_test, :);
        Y_test = Y_test(valid_test);

        % Normalisation
        [X_train, mu_X, sigma_X] = zscore(X_train);
        [Y_train, mu_Y, sigma_Y] = zscore(Y_train);

        X_test = (X_test - mu_X) ./ sigma_X;

        % Ridge regression
        beta = ridge(Y_train, X_train, lambda, 0);

        % Prediction
        Y_pred = X_test * beta(2:end);
        Y_pred = Y_pred * sigma_Y + mu_Y;
        Y_pred_all{i} = Y_pred;

        % Calculate explained variance
        res = Y_test - Y_pred;
        CVEV(i) = 1 - (var(res, 'omitnan') / var(Y_test, 'omitnan'));

        % Calculate mean squared error (MSE)
        MSE(i) = mean(res.^2, 'omitnan');
    end
end

