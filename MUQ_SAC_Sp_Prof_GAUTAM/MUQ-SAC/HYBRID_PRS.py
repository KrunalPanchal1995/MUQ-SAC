import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from scipy.optimize import minimize

def hybrid_prs_fit(X, y, degree=2, max_iter=1000, lr=0.01):
    """
    Train a Hybrid PRS model with a hyperbolic tangent term.
    
    Parameters:
    X : array-like, shape (n_samples, 1)
        Input data.
    y : array-like, shape (n_samples,)
        Target values.
    degree : int, optional
        Polynomial degree for PRS.
    max_iter : int, optional
        Maximum number of iterations.
    lr : float, optional
        Learning rate for coefficient refinement.
    
    Returns:
    dict
        Dictionary containing model coefficients and intercept.
    """
    # Polynomial Feature Transformation
    poly = PolynomialFeatures(degree,include_bias=True)
    X_poly = poly.fit_transform(X) #X_poly is a NumPy array containing the polynomial features generated from the original input features in X. It includes all polynomial combinations of the features up to the specified degree (in this case, degree=2).i.e X1,X2,X1^2, X2^2, X1X2... 
                                   #size(X_poly)-->(978,325)
    
    # Initial PRS Model
    model = LinearRegression()
    model.fit(X_poly, y) #Finding the best-fitting coefficients that minimize the error between the predicted values and the actual target values (y). y=X⋅β+ϵ, To find Beta, epsilon. If X is (978,325) and y is (978*1) then B has to be (325*1) and Epsilon (978*1). β=(X_T.X)−X_T.y
    
    # Extract Initial Coefficients
    coef = model.coef_
    #print(len(coef)) #Coeff's mean the coefficients of the X1,X2,X1^2.... not the bias(Intercept)
    intercept = model.intercept_
    print(intercept)
    
    # Define Hybrid PRS Function
    def hybrid_prs(params, X_poly, y):
        coef_tanh = params
        y_pred = X_poly @ coef_tanh + np.tanh(X_poly @ coef_tanh)
        return np.mean((y - y_pred) ** 2)
    
    # Optimize Coefficients with tanh term
    initial_params = coef
    result = minimize(hybrid_prs, initial_params, args=(X_poly, y), method='BFGS', options={'maxiter': max_iter})
    
    final_coef = result.x
    
    return {'coefficients': final_coef,'poly': poly}
def hybrid_prs_predict(X_matrix, model_params, degree=2):
    """
    Predict using the trained Hybrid PRS model.

    Parameters:
    X : array-like, shape (n_samples, 1)
        Input data.
    model_params : dict
        Dictionary containing model coefficients and intercept.
    degree : int, optional
        Polynomial degree used during training.

    Returns:
    y_pred : array-like, shape (n_samples,)
        Predicted values.
    """
    #print(X_matrix) #[[-0.57391296  0.70584836 -0.42257518 ...  0.4400418   0.83753648 0.88671626][][][]...[]] List of Lists

    y_pred = []
    for X_ in X_matrix:
        #print(X_) #A single list
        X = [X_]  #A 2D list or a 2D NumPy array.
        X_poly = model_params['poly'].transform(X)
        coef = model_params['coefficients']  #The learned coefficients (weights) from the PRS model.
        #intercept = model_params['intercept']
        y_pred.append(X_poly @ coef + np.tanh(X_poly @ coef))

    return np.asarray(y_pred).flatten()

# Example Usage

matrix_file = 'DesignMatrix.csv'
response_file= 'response_time.csv'

X = pd.read_csv(matrix_file,header=None)
y = np.log(pd.read_csv(response_file,header=None)*10000) #Normalizing the Data
#print(X)
#print(y)
#X = FeaturesTransform(X_.to_numpy())

# Ensure y is the correct shape
y = y.values.ravel() #Converts to the yData(Column) into a list [8.61497234 8.23284435 8.64134398 ... 8.61563205 8.56827424 8.57149661]
#print(y)
# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
#print(X_train) #[978 rows x 23 columns]
#print(type(X_train))
#print(X_test)  #[245 rows x 23 columns]
#print(len(y_train)) #978*1

model_params = hybrid_prs_fit(X_train, y_train, degree=2)
#print("Model Coefficients:", model_params['coefficients'])
#print("Model Intercept:", model_params['intercept'])


# Predict using the trained model
y_pred_train = hybrid_prs_predict(X_train.to_numpy(),model_params,degree=2)
y_pred_test = hybrid_prs_predict(X_test.to_numpy(), model_params, degree=2)

# Plot results
print(len(y_pred_test))
residuals_train = y_train - y_pred_train
residuals_test = y_test - y_pred_test

# Calculate maximum residual error
residuals_prs = (np.abs(y_train - y_pred_train.flatten())/y_pred_train.flatten())*100
max_residual_prs = np.max(residuals_prs)
print(f'\nMaximum Residual Error for PRS (training): {max_residual_prs}\n')

# Calculate maximum residual error
residuals_prs_test = (np.abs(y_test - y_pred_test.flatten())/y_pred_test.flatten())*100
max_residual_prs_test = np.max(residuals_prs_test)
print(f'\nMaximum Residual Error for PRS (testing): {max_residual_prs_test}\n')

# Parity plot
plt.figure(figsize=(8, 8))
# Plot the line
#plt.plot(x, y, label="y = x", color="red", linestyle="--")
plt.scatter(y_test, y_pred_test,color="black",alpha=1,)
plt.scatter(y_train, y_pred_train,color="none",alpha=1,edgecolor="green")
plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--')
plt.xlabel('Actual Values')
plt.ylabel('Predicted Values')
plt.title(f'Parity Plot\nMaximum Residual Error (training): {max_residual_prs:.2f}%\nMaximum Residual Error (testing): {max_residual_prs_test:.2f}%')
plt.show()



