import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from .duncans import DuncanResults

# Set up common plot styles for a professional look.
sns.set_style("ticks", {"grid.linestyle": "--"})
plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["axes.titlesize"] = 14
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["figure.dpi"] = 100

def _get_se_from_results(results: DuncanResults) -> float:
    """
    An internal helper function to calculate the standard error of the mean
    difference, accounting for both equal and unequal sample sizes.
    """
    ms_error = results.header["ms_error"]
    n_counts = results.means["n"].values
    
    # Check if all treatment groups have the same number of replicates.
    if len(np.unique(n_counts)) == 1:
        n = n_counts[0]
        return np.sqrt(ms_error / n)
    else:
        # If sample sizes are unequal, use the harmonic mean.
        num_treatments = len(n_counts)
        n_harm = num_treatments / np.sum(1.0 / n_counts)
        return np.sqrt(ms_error / n_harm)

def plot_bar(results: DuncanResults):
    """
    Generates a bar plot of treatment means with error bars and Duncan's grouping letters.
    
    Args:
        results (DuncanResults): A results object returned by `duncan_test`.
    """
    # A quick sanity check to ensure the input is correct.
    if not isinstance(results, DuncanResults):
        raise TypeError("Input must be a DuncanResults object from duncan_test.")

    group_data = results.grouping.copy()
    
    # Sort the data by mean for a more logical plot order, if not already sorted.
    if results.sorting != 'mean':
        group_data = group_data.sort_values("Mean", ascending=False)

    # Calculate the Standard Error (SE) using a helper function.
    se = _get_se_from_results(results)
    
    fig, ax = plt.subplots()
    
    # Create a unique color for each grouping letter for visual clarity.
    unique_groups = group_data["Group"].unique()
    palette = sns.color_palette("Set2", n_colors=len(unique_groups))
    group_to_color = {group: palette[i] for i, group in enumerate(unique_groups)}
    colors = [group_to_color[group] for group in group_data["Group"]]

    # Plot the bars with error bars.
    ax.bar(
        x=range(len(group_data)),
        height=group_data["Mean"],
        yerr=se if se and not np.isnan(se) and se > 0 else None,
        capsize=5,
        width=0.45,
        edgecolor="black",
        linewidth=0.5,
        color=colors,
        error_kw={"elinewidth": 1, "capthick": 1},
    )

    # Add the grouping letters as text labels above each bar.
    for i, (idx, row) in enumerate(group_data.iterrows()):
        text_y_pos = row["Mean"] + se + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.02
        ax.text(
            i,
            text_y_pos,
            row["Group"],
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=12,
        )

    # Finalize the plot with titles, labels, and styling.
    ax.set_title("Means with Duncan's Grouping")
    ax.set_ylabel("Mean Value")
    ax.set_xticks(range(len(group_data)))
    ax.set_xticklabels(group_data.index, fontsize=10)
    sns.despine()
    ax.yaxis.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()

def plot_cld(results: DuncanResults):
    """
    Generates a Compact Letter Display (CLD) scatter plot of means with error bars.
    
    Args:
        results (DuncanResults): A results object returned by `duncan_test`.
    """
    if not isinstance(results, DuncanResults):
        raise TypeError("Input must be a DuncanResults object from duncan_test.")

    group_data = results.grouping.copy()
    
    # Sort the data by mean in ascending order for a clean horizontal plot.
    if results.sorting == 'mean':
        group_data = group_data.sort_values("Mean", ascending=True)

    se = _get_se_from_results(results)
    
    fig, ax = plt.subplots()
    y_pos = range(len(group_data))
    
    # Use a visually pleasing color palette for the scatter points.
    colors = sns.color_palette("viridis", n_colors=len(group_data))

    # Plot the mean values as scatter points.
    ax.scatter(
        x=group_data["Mean"],
        y=y_pos,
        s=150,
        facecolor=colors,
        edgecolors="grey",
        linewidth=0.5,
        zorder=3,
    )

    # Add horizontal error bars based on the standard error.
    if se and not np.isnan(se) and se > 0:
        ax.errorbar(
            x=group_data["Mean"],
            y=y_pos,
            xerr=se,
            fmt="none",
            color="gray",
            capsize=5,
            elinewidth=1.5,
        )

    # Add the grouping letters next to each point.
    for i, (idx, row) in enumerate(group_data.iterrows()):
        text_x_pos = row["Mean"] + se + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.01
        ax.text(
            text_x_pos,
            i,
            row["Group"],
            ha="left",
            va="center",
            fontsize=12,
            fontweight="bold"
        )

    # Finalize plot details.
    ax.set_title("Compact Letter Display of Means")
    ax.set_xlabel("Mean Value")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(group_data.index, fontsize=12)
    sns.despine()
    ax.grid(axis="x", alpha=0.6)
    plt.tight_layout()

def plot_heatmap(results: DuncanResults):
    """
    Generates a heatmap of pairwise p-values from Duncan's test.
    
    Args:
        results (DuncanResults): A results object returned by `duncan_test`.
    """
    if not isinstance(results, DuncanResults):
        raise TypeError("Input must be a DuncanResults object from duncan_test.")

    # Convert the p-value matrix to a pandas DataFrame for easier plotting with Seaborn.
    pvals = results.pvals
    treatments = results.grouping.index.tolist()
    pval_df = pd.DataFrame(pvals, index=treatments, columns=treatments)

    # Create a text annotation matrix: add a '*' for significant p-values.
    annotations = pval_df.map(lambda p: f"{p:.3f}*" if p < results.header["alpha"] else f"{p:.3f}")

    fig, ax = plt.subplots(figsize=(8, 8))
    cmap = sns.diverging_palette(220, 20, as_cmap=True)

    # Generate the heatmap.
    sns.heatmap(
        pval_df,
        annot=annotations,
        fmt="", # Tells seaborn to use our string annotations.
        cmap=cmap,
        cbar_kws={"label": "P-value", "shrink": 0.75},
        vmin=0,
        vmax=1,
        ax=ax,
        square=True,
        linewidths=0.5,
        linecolor="white",
    )

    # Add a legend for the significance marker.
    legend_text = f"* indicates p < {results.header['alpha']} (significant)"
    ax.text(0.5, -0.08, legend_text, transform=ax.transAxes, ha='center', fontsize=10, style='italic')

    # Finalize plot details.
    ax.set_title("Pairwise P-values from Duncan's Test")
    plt.tight_layout()
