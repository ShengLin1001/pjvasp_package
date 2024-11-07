"""
mymetal.ml

This subpackage provides machine learning tools and utilities for materials simulations and analysis. 
It includes dataset handling, model training, data visualization, and model prediction visualization 
specifically tailored for material science applications.

Modules:
    dataset: Module for handling datasets, including functions for creating custom datasets and computing statistics.
    plot: Module for data visualization, including functions for displaying images and plots.
    model: Module for training and visualizing machine learning models.
    confusionmatrix: Module for creating confusion matrices.
"""

from mymetal.ml.dataset import (CustomDataset, get_all_data, get_mean_std)
from mymetal.ml.plot import (my_imshow)
from mymetal.ml.model import (train_model, visualize_model, visualize_model_predictions)

__all__ = ['CustomDataset', 'get_all_data', 'get_mean_std',
           'my_imshow',
           'train_model', 'visualize_model', 'visualize_model_predictions'
]
