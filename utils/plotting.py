from matplotlib import rcParams
import matplotlib.pyplot as plt
import itertools
from collections import defaultdict
import matplotlib.lines as mlines


# Constants
TEXTWIDTH_IN = 6.759374567593746
LINE_WIDTH_IN = 0.5 * TEXTWIDTH_IN


def set_params():
    """
    Set font sizes for all subsequently created plots.
    
    Args:
        None
        
    Returns:
        None
    """
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["DejaVu Sans"]
    rcParams["font.size"] = 20
    rcParams["axes.labelsize"] = 22
    rcParams["axes.titlesize"] = 26
    rcParams["legend.fontsize"] = 18
    rcParams["legend.title_fontsize"] = 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    rcParams["figure.titlesize"] = 28


def make_param_label(row):
    """
    Create a parameter label for plotting based on the method and its parameters.

    Args:
        Row-like object containing at least a method field and 
        the relevant parameter field for that method.
    
    Returns:
        Label describing the method and selected parameter value.
    """
    if row.method == "scvi":
        return f"scVI n_latent={int(row.n_latent)}"
    elif row.method == "pca":
        return f"PCA n_comp={int(row.n_components)}"
    else:
        return f"{row.method.upper()} param=unknown"


def extract_sort_key(label):
    """
    Extract method name and numeric parameter value from a plot label.

    Args:
        label: Label string such as "scVI n_latent=10" or "PCA n_comp=20".

    Returns:
        tuple: Method name and parsed numeric parameter value.
    """
    parts = label.split()
    method = parts[0]
    param_str = parts[1] if len(parts) > 1 else ""
    try:
        param_val = float(param_str.split("=")[1])
    except Exception:
        param_val = float("inf")
    return method, param_val


def plot_summary(summary, save=None, measure="V-measure"):
    """
    Plot summary statistics from a DataFrame.

    Args:
        summary :
            DataFrame containing summary statistics with columns:
            - 'method'
            - 'dropout'
            - 'median'
            - 'q10'
            - 'q90'
            - 'param_label'

        save :
            If provided, saves the plot to the specified file path.
            If None, displays the plot interactively.
        
        measure: Name of the metric shown on the y-axis.

    Returns:
        None
    """
    summary["param_label"] = summary.apply(make_param_label, axis=1)

    # Styles
    line_styles = {"scvi": "-", "pca": "--"}
    marker_cycle = itertools.cycle(["o", "s", "D", "^", "v", "P", "X", "*", "h", "8"])
    marker_cycle_map = {}

    # Labels, colors, and handles
    legend_handles = {}
    all_param_labels = sorted(summary["param_label"].unique())

    for label in all_param_labels:
        marker_cycle_map[label] = next(marker_cycle)

    # --- PLOTTING ---
    fig, ax = plt.subplots(1, 1, figsize=(TEXTWIDTH_IN, 5.0))

    for (method, param_label), group in summary.groupby(["method", "param_label"]):
        group_sorted = group.sort_values("dropout")
        ls = line_styles[method]
        marker = marker_cycle_map[param_label]

        (line,) = ax.plot(
            group_sorted["dropout"],
            group_sorted["median"],
            label=param_label,
            linestyle=ls,
            marker=marker,
            markersize=4,
        )

        ax.fill_between(
            group_sorted["dropout"],
            group_sorted["q10"],
            group_sorted["q90"],
            alpha=0.2,
        )

        if param_label not in legend_handles:
            legend_handles[param_label] = line

    # --- AXES AND LABELS ---
    baseline = summary["dropout"].min()

    xmin, xmax = ax.get_xlim()
    auto_ticks = ax.get_xticks()

    # Remove ticks before baseline and ticks too close to baseline
    min_gap = 0.02

    ticks = [
        tick for tick in auto_ticks
        if tick > baseline + min_gap and xmin <= tick <= 1.0
    ]

    ticks = [baseline] + ticks

    tick_labels = [
        "Baseline" if x == baseline else f"{x:.2f}"
        for x in ticks
    ]

    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, rotation=90)
    ax.set_xlim(left=xmin, right=1.0)
    ax.set_xlabel(r"Fraction of 0s")
    ax.set_ylabel(measure)
    #ax.set_xlim(0.85, 0.98) # zoom in
    ax.grid(True)

    # --- LEGEND ---

    # Group legend entries by method
    grouped = defaultdict(list)
    for label, handle in legend_handles.items():
        method, param_val = extract_sort_key(label)
        grouped[method].append((param_val, label, handle))

    # Compose sorted handles and labels for legend
    sorted_handles = []
    sorted_labels = []

    for method in sorted(grouped.keys()):
        # Group header line (invisible)
        dummy_line = mlines.Line2D([], [], color="none", label=method)
        sorted_handles.append(dummy_line)
        sorted_labels.append(method)

        # Sorted entries by param value
        for _, label, handle in sorted(grouped[method], key=lambda x: x[0]):
            raw_param = label.split(" ", 1)[1]  # e.g. "n_lat=10"
            name, value = raw_param.split("=")
            subscript = name.split("_")[1]
            param_only = rf"$n_{{\text{{{subscript}}}}} = {value}$"
            sorted_handles.append(handle)
            sorted_labels.append(param_only)

    # Draw legend outside the plot area
    #ax.legend(
    #    sorted_handles,
    #    sorted_labels,
    #    loc="center left",
    #    bbox_to_anchor=(1.02, 0.5),
    #    title="Method and Parameter",
    #    handlelength=2,
    #    frameon=True,
    #)

    plt.tight_layout()
    if save:
        plt.savefig(save, bbox_inches="tight", dpi=300)
    else:
        plt.show()
