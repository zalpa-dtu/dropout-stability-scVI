import re
import os
import ast


def find_pickle_files(directory):
    """
    Find all pickle result files in a directory and its subdirectories.

    Args:
        directory: Path to the directory to search.

    Returns:
        Full paths to all files ending with ".pkl".
    """
    pickle_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".pkl"):
                pickle_files.append(os.path.join(root, file))
    return pickle_files


def parse_params(param_str_or_dict):
    """
    Parse pipeline parameters from a dictionary or string representation of a dictionary.

    Args:
        param_str_or_dict: Parameter dictionary, or a string that can be
            parsed as a Python literal dictionary.

    Returns:
        Parsed parameter dictionary.
    """
    if isinstance(param_str_or_dict, dict):
        return param_str_or_dict
    try:
        return ast.literal_eval(param_str_or_dict)
    except Exception:
        return {}


def extract_params(filename):
    """
    Extract metadata and experiment parameters from a filename that follows specific patterns.

    Args:
        filename: Filename string to parse.

    Returns:
        Dictionary with extracted parameters, or None if no match.
    """

    basename = os.path.basename(filename)

    pattern_observed_counts = (
        r"^(?P<method>\w+)_(?P<simulator>\w+)_observed_counts_"
        r"(?:(?P<datatype>[^_]+)_)?"
        r"(?:(?P<dataname>[^_]+)_)?"
        r"(?P<genes>\d+)genes_"
        r"(?:(?P<factor>\d+(?:\.\d+)?(?:e[+-]?\d+)?)library_)?"
        r"(?:(?P<rate>\d+(?:\.\d+)?(?:e[+-]?\d+)?)depth_)?"
        r"(?:(?P<rate2>[\dp.]+)rate_)?"
        r"(?:(?P<mid>\d+(?:p\d+)?)mid_)?"
        r"(?:(?P<shape>m?\d+(?:p\d+)?)shape_)?"
        r"(?:(?P<alpha0>\d+(?:p\d+)?)alpha0_)?"
        r"(?:(?P<prop>\d+(?:p\d+)?)prop_)?"
        r"(?:(?P<meanlibrary>\d+(?:p\d+)?(?:e[+-]?\d+)?)meanlibrary_)?"
        r"(?:(?P<varlibrary>\d+(?:p\d+)?(?:e[+-]?\d+)?)varlibrary_)?"
        r"(?P<pct0>[\dp.]+)pct0_"
        r"\d+seed" # ignore the 110seed
        r"(?:_\{.*?\})?" # optional params dict
        r"_\{.*?\}_(?P<seed>\d+)\.pkl$"
    )

    pattern_true_gene_expr = (
        r"^(?P<method>\w+)_(?P<simulator>\w+)_true_gene_expr_"
        r"(?:(?P<datatype>[^_]+)_)?"
        r"(?:(?P<dataname>[^_]+)_)?"
        r"(?P<genes>\d+)genes_"
        r"(?:(?P<factor>\d+(?:\.\d+)?(?:e[+-]?\d+)?)library_)?"
        r"(?P<pct0>[\dp.]+)pct0_"
        r"\d+seed"
        r"(?:_\{.*?\})?"
        r"_\{.*?\}_(?P<seed>\d+)\.pkl$"
    )

    pattern_true_counts = (
        r"^(?P<method>\w+)_(?P<simulator>\w+)_true_counts_"
        r"(?P<genes>\d+)genes_"
        r"(?P<cells>\d+)cells_"
        r"(?P<pct0>[\dp.]+)pct0_"
        r"\d+seed"
        r"(?:_\{.*?\})?"
        r"_\{.*?\}_(?P<seed>\d+)\.pkl$"
    )


    def smart_cast(x):
        """
        Convert a string value to int or float when possible.

        Args:
            x: Value to convert.

        Returns:
            Converted numeric value when possible; otherwise the
            original value. Returns None if x is None.
        """
        if x is None:
            return None
        try:
            f = float(x)
            return int(f) if f.is_integer() else f
        except ValueError:
            return x

    def parse_p_float(x):
        """
        Parse numbers that may use "p" as a decimal separator.

        For example, "0p5" is converted to 0.5.

        Args:
            x: Value to parse.

        Returns:
            Parsed float when possible; otherwise the original value.
            Returns None if x is None.
        """
        if x is None:
            return None
        try:
            return float(str(x).replace("p", "."))
        except ValueError:
            return x

  
    def parse_dropout():
        """
        Match the filename against supported dropout result filename patterns.

        Args:
            None

        Returns:
            Extracted dropout experiment metadata if a filename pattern
            matches; otherwise None.
        """
        match = re.match(pattern_observed_counts, basename)
        if not match:
            match = re.match(pattern_true_gene_expr, basename)
        if not match:
            match = re.match(pattern_true_counts, basename)
            if not match:
                print(f"No pattern matched filename: {filename}")
                return None
        pct0 = parse_p_float(match.group("pct0"))
        dropout = (pct0 / 100) if isinstance(pct0, (int, float)) else pct0
        return {
            "method": match.group("method").lower(),
            "genes": smart_cast(match.group("genes")),
            "cells": smart_cast(match.groupdict().get("cells")),
            "dropout": dropout,
            "factor": parse_p_float(match.groupdict().get("factor")),
            "rate": parse_p_float(match.groupdict().get("rate") or match.groupdict().get("rate2")),
            "datatype": match.groupdict().get("datatype"),
            "dataname": match.groupdict().get("dataname"),
            "seed": int(match.group("seed")),
        }

    result = parse_dropout()
    return result


def extend_with_attributes(obj, base_dict=None):
    """
    Merges all non-private, non-callable attributes of `obj` into `base_dict`.
    Keeps nested dictionaries as values - no flattening.

    Args:
        obj: Object whose public attributes should be extracted.
        base_dict: Existing dictionary to extend.

    Returns:
        Dictionary containing base_dict values plus non-private, non-callable
        attributes from obj.
    """
    result = dict(base_dict) if base_dict else {}

    for attr in dir(obj):
        if attr.startswith("_") or callable(getattr(obj, attr)):
            continue  # Skip private and methods

        result[attr] = getattr(obj, attr)

    return result