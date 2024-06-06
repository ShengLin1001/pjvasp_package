import os
import torch
from torch.utils.data import Dataset
from PIL import Image
import torchvision

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
