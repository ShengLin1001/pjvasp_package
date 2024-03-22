import os

def print_package_structure(package_name, root_dir):
    """
    Print the directory structure of a Python package.

    Args:
        package_name (str): The name of the package.
        root_dir (str): The root directory of the package.

    Returns:
        None
    """
    print(f"Package: {package_name}")
    for root, dirs, files in os.walk(root_dir):
        level = root.replace(root_dir, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print(f"{indent}{os.path.basename(root)}/")
        subindent = ' ' * 4 * (level + 1)
        for file in files:
            print(f"{subindent}{file}")

if __name__ == "__main__":
    package_name = input("Enter the name of your package: ")
    root_dir = input("Enter the root directory of your package: ")
    print_package_structure(package_name, root_dir)

