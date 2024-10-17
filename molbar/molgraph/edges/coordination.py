import numpy as np
from scipy import special


class Coordination:

    def __init__(self) -> None:
        pass

    def _define_rij0(self, include_all=False) -> np.ndarray:
        """
        Construct a matrix containing the sum of the covalent radii of the atoms.

        Args:
            include_all (bool, optional): If all atoms in the molecule should be considered else only visible atoms. Defaults to False.

        Returns:
            np.ndarray: Matrix containing the sum of the covalent radii of the atoms.
        """

        # Get covalent radii of all atoms
        rcov = self.return_node_data(attribute="rcov", include_all=include_all)
        n_atoms = self.return_n_atoms(include_all=include_all)
        rij0 = np.ones((n_atoms, n_atoms))
        for i, ith_rcov in enumerate(rcov):
            for j, jth_rcov in enumerate(rcov):
                if j > i:
                    # Sum of covalent radii of atoms i and j
                    rij0[i, j] = ith_rcov + jth_rcov
                    rij0[j, i] = ith_rcov + jth_rcov
        return rij0

    def _return_d4(self, rij, rij0, en) -> np.ndarray:
        """
        Detemines the coordination number from the D4 dispersion correction found in
        Stefan Grimme, Jens Antony, Stephan Ehrlich, and Helge KriegJ. Chem. Phys. 132, 154104 (2010); DOI:10.1063/1.3382344


        Args:
            drij (np.ndarray): matrix of (n_atoms, n_atoms) containing the (rij-rij0)/rij0

        Returns:
            np.ndarray: returns array of shape (n_atoms,) containing the D3 coordination number of each atom.
        """

        k0 = 7.5
        k1 = 4.10451
        k2 = 19.08857
        k3 = 2 * 11.28174**2
        covCN = np.zeros(rij.shape[0])
        rij0_scaled = 4.0 / 3.0 * rij0
        for i in range(rij.shape[0]):
            for j in range(rij.shape[1]):
                if i > j:
                    r0 = rij0_scaled[i] + rij0_scaled[j]
                    r = rij[i, j]
                    abs_den = abs((en[i] - en[j]))
                    delta_en_ij = k1 * np.exp(-((abs_den + k2) ** 2) / k3)
                    erf = 1.0 + special.erf(-k0 * (r - r0) / r0)
                    tmp = delta_en_ij / 2.0 * erf
                    covCN[i] += tmp
                    covCN[j] += tmp
        return covCN

    def _get_edges_d3(self, include_all=False) -> list:
        """
        Determine whether there is an edge between two atoms or not based on the distance between them and a threshold.

        Args:
            include_all (bool, optional): If all atoms in the molecule should be considered else only visible atoms. Defaults to False.

        Returns:
            list: List of edges.
        """
        n_atoms = self.return_n_atoms(include_all=include_all)
        # Get distance matrix of all atoms
        rij = self.return_distance_matrix(include_all=include_all)
        elements = self.return_node_data(
            attribute="atomic_numbers", include_all=include_all
        )
        rij_bohr = rij / 0.52917726
        np.fill_diagonal(rij, 1.0)
        # Get matrix containing the sum of the covalent radii of the atoms
        rij0 = self.return_node_data(attribute="rcov", include_all=include_all)
        rij0_bohr = rij0 / 0.52917726
        bonds = []
        rij_norm = np.ones_like(rij, dtype=float) * np.inf
        for i in range(n_atoms):
            for j in range(n_atoms):
                if j > i:
                    value = rij[i, j] / (4.0 / 3.0 * (rij0[i] + rij0[j]))
                    # Calculate the shift of the sum of the covalent radii of the atoms based on the D3 coordination number. More Information can be found in the cited paper.
                    rij_norm[i, j] = value
                    rij_norm[j, i] = value

        rij_norm = np.where(rij_norm < 1.0, rij_norm, np.inf)
        bonds = list(zip(*np.where(rij_norm < 1.0)))
        elements = self.return_node_data(attribute="elements", include_all=include_all)
        rij0[np.where(elements == "H")[0]] += 0.3
        return bonds, rij0

    def _get_edges_gfnff(self, include_all=False) -> list:
        """
        Determine whether there is an edge between two atoms or not based on the distance between them and a threshold.
        More information can be3 found in the SI of S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665.

        Args:
            include_all (bool, optional): If all atoms in the molecule should be considered else only visible atoms. Defaults to False.

        Returns:
            list: List of edges.
        """
        n_atoms = self.return_n_atoms(include_all=include_all)
        # Get distance matrix of all atoms
        rij = self.return_distance_matrix(include_all=include_all)
        rij_bohr = rij / 0.52917726
        np.fill_diagonal(rij, 1.0)
        # Get matrix containing the sum of the covalent radii of the atoms
        rij0 = self.return_node_data(attribute="rcov", include_all=include_all)
        rij0_bohr = self.return_node_data(attribute="rnorm", include_all=include_all)
        # rij0_bohr = rij0/0.52917726
        # Get tabulated fit parameters for each element.
        cn_fak = self.return_node_data(attribute="cn_fak", include_all=include_all)
        per_fak = self.return_node_data(attribute="per_fak", include_all=include_all)
        en = self.return_node_data(attribute="en_pauling", include_all=include_all)
        # Get the tabulated coordination number for elements.
        normcn = self.return_node_data(attribute="normcn", include_all=include_all)
        bonds = []
        for i in range(n_atoms):
            for j in range(n_atoms):
                if j > i:
                    # Calculate the shift of the sum of the covalent radii of the atoms based on the D3 coordination number. More Information can be found in the cited paper.
                    r_i = rij0_bohr[i] + cn_fak[i] * normcn[i]
                    r_j = rij0_bohr[j] + cn_fak[j] * normcn[j]
                    delta_en = np.abs(en[i] - en[j])
                    k1 = 5e-3 * (per_fak[i, 0] + per_fak[j, 0])
                    k2 = 5e-3 * (per_fak[i, 1] + per_fak[j, 1])
                    factor = 1 - k1 * delta_en - k2 * delta_en**2
                    rij_ref = (r_i + r_j) * factor
                    if rij_bohr[i, j] < 1.25 * rij_ref:
                        bonds.append((i, j))
        elements = self.return_node_data(attribute="elements", include_all=include_all)
        rij0[np.where(elements == "H")[0]] += 0.3
        return bonds, rij0
