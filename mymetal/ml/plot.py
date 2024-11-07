"""
plot module

This module provides functions for data visualization, including displaying images and plots.

Functions:
    - my_imshow: Display image for Tensor.
"""

import numpy as np
import matplotlib.pyplot as plt

def my_imshow(inp, title=None):
    """Display image for Tensor."""
    inp = inp.numpy().transpose((1, 2, 0))  # 转换为 NumPy 数组并调整维度顺序
    inp = inp.squeeze()  # 移除单通道
    mean = np.array([0.5])
    std = np.array([0.5])
    inp = std * inp + mean  # 还原归一化
    inp = np.clip(inp, 0, 1)  # 确保数值在 0 和 1 之间
    plt.imshow(inp, cmap='gray')  # 显示灰度图像
    if title is not None:
        plt.title(title)
    plt.pause(0.001)  # pause a bit so that plots are updated
    plt.show()  # 显示图像
