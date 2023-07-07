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
    Modified version of the sklearn one.
    It is also known as the Variance Ratio Criterion.

    The score is defined as ratio of the sum of between-cluster dispersion and
    of within-cluster dispersion.

    Read more in the :ref:`User Guide <calinski_harabasz_index>`.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        A list of ``n_features``-dimensional data points. Each row corresponds
        to a single data point.

    labels : array-like of shape (n_samples,)
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
    matrix_mean = np.mean(np.mean(X))
    outliers_penalty = matrix_mean**2 * 9

    check_number_of_labels(n_labels, n_samples)
    extra_disp, intra_disp = 0.0, 0.0
    mean = np.mean(X, axis=0)
    for k in range(n_labels):
        
        cluster_k = X[labels == k]
        mean_k = np.mean(cluster_k, axis=0)

        extra_disp += len(cluster_k) * np.sum((mean_k - mean) ** 2)
        if len(cluster_k) <= 1: # outliers
            intra_disp += outliers_penalty # instead of 0
        else: 
            intra_disp += np.sum((cluster_k - mean_k) ** 2)

    return (
        1.0
        if intra_disp == 0.0
        else extra_disp * (n_samples - n_labels) / (intra_disp * (n_labels - 1.0))
    )

    
