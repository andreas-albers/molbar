import numpy as np
from scipy.spatial import distance_matrix
from scipy.stats import rankdata


class RelativeConfIndex:

    def get_relconf_index(
        self, geometry: np.ndarray, charges: np.ndindex, integer=True
    ) -> int or float:
        """
        Calculates the relative configuration index for a given geometry.

        Args:
            geometry (np.ndarray): Geometry of the molecule.
            charges (np.ndindex):  Charges of the atoms in the molecule.
            integer (bool, optional): Whether to return the index as an integer. Defaults to True.

        Returns:
            int or float: Relative configuration index.
        """
        n_atoms = geometry.shape[0]

        if n_atoms == 1:

            return 0.0

        else:

            rij = distance_matrix(geometry, geometry)

            np.fill_diagonal(rij, 1)

            p2 = np.outer(charges, charges.T)

            qij = np.multiply(p2 * 100, 1 / rij)

            B = qij.sum()

        if integer == True:

            return int(np.round(B, 3))

        else:

            return B
