import matplotlib.pyplot as plt
import numpy as np

# Define a mapping for interaction names
NAMING_MAP = {
    "HBAcceptor": "HBA",
    "HBDonor": "HBD",
    "ImplicitHBAcceptor": "HBA",
    "ImplicitHBDonor": "HBD",
}


def get_interactions(df):
    "Using interaction df to get a set of interactions."

    interaction_set = set()
    interaction_set.update(
        (
            each_interaction[0],
            each_interaction[1],
            NAMING_MAP.get(each_interaction[-1], each_interaction[-1]),
        )
        for each_interaction in set(df[0].index)
    )

    return interaction_set


def confusion_matrix(explicit_set, implicit_set):
    """
    Plot a confusion matrix of explicit and implicit interactions.
    """
    intersection = explicit_set.intersection(implicit_set)
    only_explicit = explicit_set - intersection
    only_implicit = implicit_set - intersection

    return np.array([[len(intersection), len(only_explicit)], [len(only_implicit), 0]])


def plot_confusion_matrix(
    matrix,
    title="Confusion Matrix of Explicit and Implicit Hbond Interactions",
):
    fig, ax = plt.subplots(dpi=300)
    cax = ax.imshow(matrix, cmap="Oranges")
    fig.colorbar(cax)

    ax.set_xticks(np.arange(2))
    ax.set_yticks(np.arange(2))
    ax.set_xticklabels(["Detected by", "Not Detected by"])
    ax.set_yticklabels(["Detected by", "Not Detected by"])

    for (i, j), val in np.ndenumerate(matrix):
        ax.text(j, i, val, ha="center", va="center")

    ax.set_title(title)
    ax.set_xlabel("Implicit Method")
    ax.set_ylabel("Explicit Method")

    return fig, ax


def tanimoto_coefficient(set1, set2):
    """
    Compute the Tanimoto coefficient between two sets.
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    if union == 0:
        return 1  # No interactions by either method (correct)
    return intersection / union


def positive_predictive_value(true_set, predicted_set):
    """
    Compute the Positive Predictive Value (PPV) between two sets.
    """
    true_positives = len(true_set.intersection(predicted_set))
    false_positives = len(predicted_set - true_set)
    if true_positives + false_positives == 0:
        return 1  # No positive predictions (correct)
    return true_positives / (true_positives + false_positives)


def sensitivity(true_set, predicted_set):
    """
    Compute the Sensitivity between two sets.
    """
    true_positives = len(true_set.intersection(predicted_set))
    false_negatives = len(true_set - predicted_set)
    if true_positives + false_negatives == 0:
        return 1  # No actual positives (correct)
    return true_positives / (true_positives + false_negatives)


def tanimoto_coefficient_by_confusion_matrix(matrix):
    """
    Compute the Tanimoto coefficient from a confusion matrix.
    """
    intersection = matrix[0, 0]
    union = np.sum(matrix) - matrix[1, 1]
    if union == 0:
        return 1  # No interactions by either method (correct)
    return intersection / union
