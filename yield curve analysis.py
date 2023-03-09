# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 12:34:13 2023

@author: Gilberto
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize

# Define Vasicek factor processes
def dx(t, x, k1, v1, dw1):
    return k1 * (v1 - x) * t + np.sqrt(v1) * dw1

def dy(t, y, k2, v2, dw2):
    return k2 * (v2 - y) * t + np.sqrt(v2) * dw2

# Define zero-coupon bond price function
def zcb_price(t, T, k1, v1, k2, v2, x, y):
    B1 = (1 - np.exp(-k1*(T-t)))/k1
    B2 = (1 - np.exp(-k2*(T-t)))/k2
    A1 = np.exp(((v1/(2*k1**2)) * (B1 - T + t)**2) - ((x**2)/(2*v1)) * B1)
    A2 = np.exp(((v2/(2*k2**2)) * (B2 - T + t)**2) - ((y**2)/(2*v2)) * B2)
    return A1 * A2 * np.exp(-x*B1 - y*B2)

# Define model-implied swap rate function
def swap_rate(t, swap_tenors, k1, v1, k2, v2, x, y):
    p0t = zcb_price(0, t, k1, v1, k2, v2, x, y)
    pn = np.sum([zcb_price(t, swap_tenors[i], k1, v1, k2, v2, x, y) for i in range(len(swap_tenors))])
    return (p0t - pn)/np.sum([zcb_price(0, swap_tenors[i], k1, v1, k2, v2, x, y) for i in range(len(swap_tenors))])

# Define objective function for calibration
def objective_function(params, t, swap_tenors, swap_rates):
    k1, v1, k2, v2, x, y = params
    error = 0
    for i in range(len(t)):
        p1 = swap_rate(t[i], swap_tenors, k1, v1, k2, v2, x[i], y[i])
        p10 = swap_rate(t[i] + 10, swap_tenors, k1, v1, k2, v2, x[i], y[i])
        error += (p1 - swap_rates[i][0])**2 + (p10 - swap_rates[i][1])**2
    return error

# Define rolling window for calibration
window_size = 5
step_size = 1

# Define Vasicek parameters from Duarte et al. (2006)
k1 = 0.0009503
v1 = 0.0548290
k2 = 0.0113727
v2 = 0.4628664

# Define initial values for x and y
x0 = np.zeros(window_size)
y0 = np.zeros(window_size)

# Load swap curve data
swap_curve = pd.read_csv("swap_curve.csv", index_col=0)
swap_tenors = swap_curve.index.values
swap_rates = swap_curve.values

# Perform rolling window calibration
for i in range(window_size, len(swap_curve), step_size):
    # Getthe calibration window
    calibration_window = swap_curve[i-window_size:i]
# Define the target swap rates for 1-year and 10-year tenors
target_swap_rates = np.array([calibration_window.loc[1]['Rate'], calibration_window.loc[10]['Rate']])

# Define the objective function to minimize (difference between model implied swap rates and target swap rates)
objective_function = lambda x: np.linalg.norm(model_implied_swap_rates(calibration_window.index, x[0], x[1], x[2], x[3], x[4], x[5]) - target_swap_rates)

# Minimize the objective function to find the optimal parameters
initial_guess = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01] # starting point for optimization
optimal_params = minimize(objective_function, initial_guess).x

# Store the optimal parameters for this calibration window
optimal_params_df.loc[calibration_window.index[-1]] = optimal_params


