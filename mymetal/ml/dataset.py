"""
dataset module

This module provides functions for handling datasets, including creating custom datasets and computing statistics.

Functions:
    - CustomDataset: Custom dataset class for handling image data.
    - get_all_data: Get all image file paths and their corresponding labels.
    - get_mean_std: Compute mean and standard deviation for training data.
"""

import os
import torch
from torch.utils.data import Dataset
from PIL import Image
import torchvision
from torchvision.datasets import ImageFolder

class CustomDataset(Dataset):
    def __init__(self, img_paths: list=None, labels: list=None, transform: torchvision.transforms.transforms.Compose=None,
                  color: str="L", categories: list = None):
        """
        img_paths = ["./...", "./..."]
        labels = ['','']
        data_transforms = {
        'train': transforms.Compose([
            transforms.Resize((224, 224)),  # 直接调整图像大小
            transforms.RandomHorizontalFlip(),
            transforms.ToTensor(),
            transforms.Normalize([0.5], [0.5])  # 归一化灰度图像
        ]),
        'val': transforms.Compose([
            transforms.Resize((224, 224)),  # 直接调整图像大小
            transforms.ToTensor(),
            transforms.Normalize([0.5], [0.5])  # 归一化灰度图像
        ]),
        'predict': transforms.Compose([
            transforms.Resize((224, 224)),  # 直接调整图像大小
            transforms.ToTensor(),
            transforms.Normalize([0.5], [0.5])  # 归一化灰度图像
        ]),
        }
        """
        self.img_paths = img_paths
        self.labels = labels
        self.transform = transform
        self.label_to_idx = {label: idx for idx, label in enumerate(categories)}
        self.color = color

    def __len__(self):
        '''
        len(customdataset)
        '''
        return len(self.img_paths)

    def __getitem__(self, idx):
        '''
        customdataset[0]
        '''
        img_path = self.img_paths[idx]
        image = Image.open(img_path).convert(self.color)  # 转换为灰度图像/彩色
        label = self.labels[idx]
        label_idx = self.label_to_idx[label]

        if self.transform:
            image = self.transform(image)

        return image, img_path, label_idx,  label
    

def get_all_data(data_dir: str=None, categories: list=None, file_end: tuple=('.png', '.jpg', '.jpeg', '.bmp', '.gif') ) -> tuple:
    '''
    Get all image file paths and their corresponding labels from the specified directory structure.

    The function traverses the given directory for each category and collects paths of all images
    with specified file extensions along with their labels.

    Parameters
    ----------
    data_dir : str, optional
        The root directory containing category subdirectories (default is None).
    categories : list, optional
        List of category names corresponding to subdirectory names in data_dir (default is None).

    Returns
    -------
    tuple
        A tuple containing two lists:
            - all_data: list of str
                List of image file paths.
            - all_labels: list of str
                List of labels corresponding to the categories of the images.

    Usage
    -----
    data_dir = './data'
    categories = ['CT_COVID', 'CT_NonCOVID']
    all_data, all_labels = get_all_data(data_dir, categories)
    
    Example
    -------
    >>> data_dir = './data'
    >>> categories = ['CT_COVID', 'CT_NonCOVID']
    >>> all_data, all_labels = get_all_data(data_dir, categories)
    >>> print(all_data)
    ['./data/CT_COVID/image1.png', './data/CT_COVID/image2.jpg', ...]
    >>> print(all_labels)
    ['CT_COVID', 'CT_COVID', 'CT_NonCOVID', ...]
    '''
    all_data = []
    all_labels = []

    for category in categories:
        category_path = os.path.join(data_dir, category)
        for img_name in os.listdir(category_path):
            img_path = os.path.join(category_path, img_name)
            if img_path.lower().endswith(file_end):
                all_data.append(img_path)
                all_labels.append(category)
    return all_data, all_labels


def get_mean_std(train_data: ImageFolder = None, colors: str = 'RGB') -> tuple:
    '''
    Compute mean and standard deviation for training data.

    Parameters:
    -----------
    train_data : torchvision.datasets.ImageFolder
        Custom Dataset or ImageFolder containing the training data.
    colors : str, optional
        Color mode of the images ('RGB' or 'gray'). Default is 'RGB'.

    Returns:
    --------
    tuple
        A tuple containing the mean and standard deviation of the training data.
        (mean, std), where mean and std are lists of floats.

    Example:
    --------
    >>> data_dir = './data_prepared/train'
    >>> train_dataset = ImageFolder(data_dir, transform=mytransform)
    >>> mean, std = get_mean_std(train_dataset, colors='gray')
    >>> print(mean, std)
    '''
    
    # Input checks
    if train_data is None:
        raise ValueError("train_data cannot be None. Please provide a valid ImageFolder dataset.")
    
    if not isinstance(train_data, ImageFolder):
        raise TypeError("train_data must be an instance of torchvision.datasets.ImageFolder.")
    
    if colors not in ['RGB', 'rgb', 'g', 'gray', 'G', 'Gray']:
        raise ValueError("Invalid value for colors. Choose from 'RGB', 'rgb', 'g', 'gray', 'G', 'Gray'.")

    print('Compute mean and variance for training data.')
    print(len(train_data))
    
    if colors in ['RGB', 'rgb']:
        temp = 3
    elif colors in ['g', 'gray', 'G', 'Gray']:
        temp = 1
    
    train_loader = torch.utils.data.DataLoader(
        train_data, batch_size=1, shuffle=False, num_workers=0, pin_memory=True)
    
    mean = torch.zeros(temp)
    std = torch.zeros(temp)
    
    for X, _ in train_loader:
        for d in range(temp):
            mean[d] += X[:, d, :, :].mean()
            std[d] += X[:, d, :, :].std()
    
    mean.div_(len(train_data))
    std.div_(len(train_data))
    
    return (list(mean.numpy()), list(std.numpy()))
