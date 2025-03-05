function cvev = cal_cvev(Y_true, Y_pred)
    if var(Y_true) == 0
        cvev = NaN;
    else
        cvev = 1 - var(Y_true - Y_pred) / var(Y_true);
    end
end
