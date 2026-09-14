import os
import pandas as pd
import scanpy as sc


def find_datasets(dir, contains=None, not_contains=None):
    """
    Find all dataset files in the given directory.

    Args:
        dir (str): The directory to search for dataset files.
        contains list(str): A list of strings that the dataset file names should contain.

    Returns:
        list: A list of paths to the dataset files fulfilling the criteria.
    """
    dataset_files = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if (
                file.endswith(".csv")
                or file.endswith(".json")
                or file.endswith(".h5ad")
            ):
                if contains is None or any(c in file for c in contains):
                    if not_contains is None or all(c not in file for c in not_contains):
                        dataset_files.append(os.path.join(root, file))
    return dataset_files


def find_pickle_files(directory):
    """
    Find all pickle files in a directory and its subdirectories.

    Args:
        directory: Path to the directory to search.

    Returns:
        list: List of full paths to files ending with ".pkl".
    """
    pickle_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".pkl"):
                pickle_files.append(os.path.join(root, file))
    return pickle_files


def load_dataset(file_path):
    """
    Load a dataset from the given .h5ad file.

    Args:
        file_path: The path to the .h5ad dataset file.

    Returns:
        Loaded AnnData object.
    """
    if not file_path.endswith(".h5ad"):
        raise ValueError(f"Expected an .h5ad file, got: {file_path}")

    return sc.read_h5ad(file_path)
