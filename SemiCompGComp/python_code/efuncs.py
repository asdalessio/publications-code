import numpy as np
from delicatessen.estimating_equations import ee_mlogit


def predict_from_mlogit(theta, X, n_y_vals):
    n_x_vals = X.shape[1]    # Number of columns / predictors in X
    denom = 1                # Default value for denominator
    exp_pred_y = []          # Storage for the expected values of Y
    start_index = 0          #

    for i in range(1, n_y_vals):                      # Looping over all columns of Y
        end_index = start_index + n_x_vals            # ... get the current end_index
        beta_i = theta[start_index: end_index]        # ... grab the corresponding beta's for Y column
        pred_y = np.exp(np.dot(X, beta_i))            # ... generate predicted value of Y column
        exp_pred_y.append(pred_y)                     # ... store the particular predicted values for Y
        denom = denom + pred_y                        # ... update the denominator with summation
        start_index = end_index                       # ... update start_index to current end_index

    ygt0_hats = np.asarray(exp_pred_y) / denom
    y0_hat = 1 - np.sum(ygt0_hats, axis=0)
    return np.vstack([y0_hat, ygt0_hats])


def ef_multistate_ice_gcomputation(theta, y_array, X_array, Xa_array, c_array):
    k_times = len(y_array)
    y_vals = 3
    mu_y = theta[:y_vals]
    beta_s = theta[y_vals:]

    # Initial time model
    x_str = 0
    x_end = X_array[0].shape[1]*2
    yhat = y_array[0]

    ee_fms = []
    for k in range(k_times):
        beta_k = beta_s[x_str: x_end]

        # Fitting nuisance model
        ee_fmk = ee_mlogit(beta_k, X=X_array[k], y=yhat) * (1-c_array[k])
        ee_fms.append(ee_fmk)

        # Predictions from the nuisance model
        yhat = predict_from_mlogit(theta=beta_k, X=Xa_array[k], n_y_vals=y_vals)
        yhat = yhat.T
        if k != k_times-1:
            y_prev = y_array[k + 1]
            yhat[y_prev[:, -1] == 1] = [0, 0, 1]
            # Updating indices
            x_str = x_end
            x_end = x_str + X_array[k].shape[1]*2

    ee_muj = []
    for j in range(y_vals):
        ee_muj.append(yhat[:, j] - mu_y[j])

    # Return stacked estimating functions
    return np.vstack(ee_muj + ee_fms)
