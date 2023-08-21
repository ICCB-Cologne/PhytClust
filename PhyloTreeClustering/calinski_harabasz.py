import functools
import numpy as np

from sklearn.utils import check_X_y
from sklearn.preprocessing import LabelEncoder


def check_number_of_labels(n_labels, n_samples):
        """Check that number of labels are valid.

        Parameters
        ----------
        n_labels : int
            Number of labels.

        n_samples : int
            Number of samples.
        """
        if not 1 < n_labels < n_samples:
            raise ValueError(
                "Number of labels is %d. Valid values are 2 to n_samples - 1 (inclusive)"
                % n_labels
            )

def calinski_harabasz_score(X, labels):
    """Compute the Calinski and Harabasz score.

    Parameters
    ----------
    X : distance matrix of shape (n_samples, n_samples)
        Each row/column represents a node of the phylogenetic tree

    labels : array-like of shape (n_samples, labels)
        Predicted labels for each sample.

    Returns
    -------
    score : float
        The resulting Calinski-Harabasz score.

    References
    ----------
    .. [1] `T. Calinski and J. Harabasz, 1974. "A dendrite method for cluster
       analysis". Communications in Statistics
       <https://www.tandfonline.com/doi/abs/10.1080/03610927408827101>`_
    """

    X, labels = check_X_y(X, labels)
    le = LabelEncoder()
    labels = le.fit_transform(labels)

    n_samples, _ = X.shape
    n_labels = len(le.classes_)
    matrix_mean = np.mean(X)/2
    outliers = X[:, labels <= 1]
    outliers = outliers[labels <= 1]
    outlier_penalty = (matrix_mean**2)*(len(outliers)/len(X)) #For all clusters with just one terminal node

    check_number_of_labels(n_labels, n_samples)
    within_cluster, between_cluster = 0.0, 0.0
    overall_ss = (np.sum(X**2))/(n_samples * 2)

    for k in range(n_labels):
        
        cluster_k = X[:, labels == k]
        cluster_k = cluster_k[labels == k]
        if len(cluster_k) <= 1:
            within_cluster += outlier_penalty 
        else: 
            within_cluster += (np.sum(cluster_k**2))/(len(cluster_k) * 2)

    return (
        1.0
        between_cluster * (n_samples - n_labels) / (within_cluster * (n_labels - 1.0))
    )