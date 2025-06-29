import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2


def multiple_histograms(data_lists, file_path, titles, xlabel, ylabel, bins=20, color='skyblue', edgecolor='black'):
    """
    Saves multiple histograms (one per list of data) in subplots, all with the same X and Y axes.

    :param data_lists: List of lists, each containing data for one histogram.
    :param file_path: Path to save the image with the subplots.
    :param titles: List of titles for each subplot.
    :param xlabel: Label for the x-axis.
    :param ylabel: Label for the y-axis.
    :param bins: Number of bins for the histograms.
    :param color: Color of the histogram bars.
    :param edgecolor: Color of the edges of the histogram bars.
    """
    # Number of subplots
    num_plots = len(data_lists)

    # Create subplots, adjusting the layout to fit all subplots in one figure
    fig, axes = plt.subplots(1, num_plots, figsize=(8 * num_plots, 6), sharex=True, sharey=True)

    # If there's only one plot, axes is not a list, so convert it to a list
    if num_plots == 1:
        axes = [axes]

    # Iterate over the data lists and titles to plot each histogram
    for i, (data, title) in enumerate(zip(data_lists, titles)):
        sns.histplot(data, bins=bins, ax=axes[i], kde=True)

        axes[i].set_title(f"{title} Gene Lengths")
        axes[i].set_xlabel(xlabel)
        axes[i].set_ylabel(ylabel)
        axes[i].grid(axis='y', linestyle='--', alpha=0.7)

    # Adjust layout so labels do not overlap
    plt.tight_layout()

    # Save the figure with all the subplots
    plt.savefig(file_path)
    print(f"Multiple histograms saved at '{file_path}'")
    plt.close()


def histogram(data, file_path, title, xlabel, ylabel, bins=20, color='skyblue', edgecolor='black'):
    """
    Saves a histogram of the provided data to the specified file path.

    :param data: List of data values for the histogram.
    :param file_path: Path to save the histogram image.
    :param title: Title of the histogram.
    :param xlabel: Label for the x-axis.
    :param ylabel: Label for the y-axis.
    :param bins: Number of bins for the histogram.
    :param color: Color of the histogram bars.
    :param edgecolor: Color of the edges of the histogram bars.
    """
    plt.figure(figsize=(8, 6))
    sns.histplot(data, bins=bins, kde=True)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(file_path)
    print(f"Histogram saved at '{file_path}'")
    plt.close()


def no_line_histogram(data, file_path, title, xlabel, ylabel, bins=20, color='skyblue', edgecolor='black'):
    """
    Saves a histogram of the provided data to the specified file path.

    :param data: List of data values for the histogram.
    :param file_path: Path to save the histogram image.
    :param title: Title of the histogram.
    :param xlabel: Label for the x-axis.
    :param ylabel: Label for the y-axis.
    :param bins: Number of bins for the histogram.
    :param color: Color of the histogram bars.
    :param edgecolor: Color of the edges of the histogram bars.
    """
    plt.figure(figsize=(8, 6))
    plt.hist(data, bins, color='skyblue', edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(file_path)
    print(f"Histogram saved at '{file_path}'")
    plt.close()


def venn_diagram(sets, set_labels, file_path, title):
    """
    Saves a Venn diagram comparing two gene sets.

    :param sets: List of two sets to be compared.
    :param set_labels: List of two labels for the sets.
    :param file_path: Path to save the Venn diagram image (default: 'venn_diagram.png').
    """
    # Check that exactly two sets are provided
    if len(sets) != 2:
        raise ValueError("Only 2 sets are supported for the Venn diagram.")

    # Check that the number of sets matches the number of labels
    if len(sets) != len(set_labels):
        raise ValueError("The number of sets must match the number of set labels.")

    # Create the Venn diagram for two sets
    plt.figure(figsize=(8, 8))
    venn = venn2(sets, set_labels)

    # Customize colors
    venn.get_patch_by_id('10').set_color('#ff9999')
    venn.get_patch_by_id('01').set_color('#66b3ff')
    venn.get_patch_by_id('11').set_color('#99ff99')

    # Set title and grid
    plt.title(title, fontsize=14)
    plt.grid(axis='both', linestyle='--', linewidth=0.5, alpha=0.7)

    # Save the Venn diagram image
    plt.tight_layout()
    plt.savefig(file_path, dpi=300)
    print(f"Venn diagram saved at '{file_path}'")
    plt.close()


def bar_plot(categories, values, file_path, title, xlabel, ylabel):
    """
    Saves a bar plot comparing different categories.

    :param categories: List of categories for the x-axis.
    :param values: List of values corresponding to each category.
    :param file_path: Path to save the bar plot image.
    :param title: Title of the bar plot.
    :param xlabel: Label for the x-axis.
    :param ylabel: Label for the y-axis.
    """
    # Check that the number of categories matches the number of values
    if len(categories) != len(values):
        raise ValueError("The number of categories must match the number of values.")

    # Create bar plot
    plt.figure(figsize=(8, 6))

    plt.bar(categories, values, color=['#ff9999', '#66b3ff', '#99ff99'])  # Default colors

    # Set the plot title and labels
    plt.title(title, fontsize=14)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Set grid and layout
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    # Save the plot
    plt.savefig(file_path, dpi=300)
    print(f"Bar plot saved at '{file_path}'")
    plt.close()


def pie_chart(categories, values, file_path, title):
    """
    Saves a pie chart showing the distribution of different categories.

    :param categories: List of category labels for the pie chart.
    :param values: List of values corresponding to each category.
    :param file_path: Path to save the pie chart image.
    :param title: Title of the pie chart.
    """
    # Check that the number of categories matches the number of values
    if len(categories) != len(values):
        raise ValueError("The number of categories must match the number of values.")

    colors = ['#ff9999', '#66b3ff', '#99ff99']

    # Create pie chart
    plt.figure(figsize=(8, 6))
    plt.pie(values, labels=categories, colors=colors, autopct='%1.1f%%', startangle=140)

    # Set title
    plt.title(title, fontsize=14)

    # Save the pie chart
    plt.tight_layout()
    plt.savefig(file_path, dpi=300)
    print(f"Pie chart saved at '{file_path}'")
    plt.close()
