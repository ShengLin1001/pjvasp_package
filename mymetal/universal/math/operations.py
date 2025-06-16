import numpy as np
from scipy.integrate import cumulative_trapezoid

def get_integration(x_values: list=None, y_values: list=None, x_target: float=None) -> float:
    """
    Integrate y-values with respect to x-values using the trapezoidal rule up to a specified x_target.
    
    Args:
        x_values (list or np.array): Sequence of x values.
        y_values (list or np.array): Corresponding sequence of y values.
        x_target (float): The x value up to which the integration is to be performed.
        
    Returns:
        float: The integral value up to x_target.
    """
    
    x_values = np.array(x_values)
    y_values = np.array(y_values)
    
    sorted_indices = np.argsort(x_values)
    x_values = x_values[sorted_indices]
    y_values = y_values[sorted_indices]

    if x_target < x_values[0] or x_target > x_values[-1]:
        raise ValueError(f"x_target {x_target} is out of the range of x_values ({x_values[0]}, {x_values[-1]})")

    idx = np.searchsorted(x_values, x_target)

    integral_value = 0.0
    for i in range(1, idx):
        integral_value += (x_values[i] - x_values[i - 1]) * (y_values[i] + y_values[i - 1]) / 2
    
    if idx > 0:
        integral_value += (x_target - x_values[idx - 1]) * (y_values[idx - 1] + (y_values[idx] if idx < len(x_values) else y_values[idx - 1])) / 2

    return integral_value


def get_integration_curve(x_values: list = None, y_values: list = None, dx=1.0, axis=-1, initial=0) -> list:
    """
    Compute the cumulative integral of y-values with respect to x-values using the trapezoidal rule.
    
    Args:
        x_values (list or np.array): Sequence of x values.
        y_values (list or np.array): Corresponding sequence of y values.
        
    Returns:
        list: The cumulative integral curve.
    """
    
    x_values = np.array(x_values)
    y_values = np.array(y_values)

    sorted_indices = np.argsort(x_values)
    x_values = x_values[sorted_indices]
    y_values = y_values[sorted_indices]

    integral_curve = cumulative_trapezoid(y_values, x=x_values, dx=dx, axis=axis, initial=initial)
    
    return integral_curve
