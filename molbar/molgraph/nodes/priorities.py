import numpy as np
from scipy.stats import rankdata
from molbar.helper.rank import Rank


class Priorities:

    def cluster_priorities(self, cluster):
        n_cluster = len(cluster)
        all_priorities = []
        for nodes in cluster:
            self.set_nodes_visible(visible_nodes=nodes)
            all_priorities.append(list(self.return_node_data(attribute="priorities")))
        # Find the maximum number of atoms in any fragment
        max_n_nodes = 0
        for j in range(n_cluster):
            max_n_nodes = max(max_n_nodes, len(all_priorities[j]))
        # Create a list of filler values to pad the fragment priorities to the same length
        filler = [0] * max_n_nodes
        # Pad the priorities of atoms in each fragment and calculate the fragment priorities
        same_length_priorities = [
            nodes_priorities[:max_n_nodes] + filler[len(nodes_priorities) :]
            for nodes_priorities in all_priorities
        ]
        return Rank(same_length_priorities).priorities

    def get_densities(self, eps, C, thrs=8):
        occ_vec = np.where(eps > 1e-3, 2, 0)
        P = C @ np.diag(occ_vec) @ C.T
        densities = np.diag(P)
        # normed_densities = densities / np.linalg.norm(densities)
        rounded_normed_densities = np.round(densities, thrs)
        return rounded_normed_densities

    def calculate_priorities(self, type="topology", thrs=8, include_all=False):
        """
        Returns the priority of each atom based on the mass spheres around it.

        Returns:
            np.ndarray: Array of priority values for each atom.
        """
        if include_all:
            self.set_all_visible()
        if type == "topology":
            spectrum = self.topology_spectrum
            orbitals = self.topology_orbitals
        elif type == "topography":
            spectrum = self.topography_spectrum
            orbitals = self.topography_orbitals
        densities = self.get_densities(spectrum, orbitals, thrs)
        priorities = rankdata(densities, method="min")
        self.add_node_data(attribute="densities", new_data=densities)
        self.add_node_data(attribute="priorities", new_data=priorities)

    def _adjust_priorities_for_pseudochirality(
        self, fragment_indices, absconf_index
    ) -> None:
        """
        Adjusts the priorities of the nodes based on the core and absolute configuration indices to account
        for pseudo chirality.

        Args:
            fragment_indices (list): indices of nodes in the fragment
            absconf_index (int): absolute configuration index of the fragment
        """
        if absconf_index == 1:
            priority_shift = 0.49
        else:
            priority_shift = -0.49
        for index in fragment_indices:
            self.nodes[index]["priorities"] += priority_shift

    def _get_atomic_spheres(self, atom: int) -> list:
        """
        NOT USED ANYMORE
        Constructs the atomic spheres around each atom based on atomic indices

        Args:
            atom (int): atomic index for which atomic spheres are to be obtained
        Returns:
            list: atomic spheres for atom, nth row -> atoms in molecule with distance n to atom.
        """

        # Get the distances of all the atoms in the molecule from the given atom.

        distances = self.distance_matrix[atom]

        # Create the atomic spheres for the given atom.
        spheres = [
            [node for node in self.visible_nodes if distances[node] == value]
            for value in np.unique(list(distances.values()))
            if value != np.inf
        ]

        # Update the atomic spheres based on the double and triple bonds.
        for i in range(len(spheres) - 1):
            for double_bond in self.double_bonds:
                if (double_bond[0] in spheres[i]) & (double_bond[1] in spheres[i + 1]):
                    spheres[i + 1].append(double_bond[1])
                if (double_bond[1] in spheres[i]) & (double_bond[0] in spheres[i + 1]):
                    spheres[i + 1].append(double_bond[0])
            for triple_bond in self.triple_bonds:
                if (triple_bond[0] in spheres[i]) & (triple_bond[1] in spheres[i + 1]):
                    spheres[i + 1] += 2 * [triple_bond[1]]
                if (triple_bond[1] in spheres[i]) & (triple_bond[0] in spheres[i + 1]):
                    spheres[i + 1] += 2 * [triple_bond[0]]
        return spheres

    def _adjust_size_mass_spheres(self, mass_spheres_for_each_atom: list) -> list:
        """
        NOT USED ANYMORE
        Fills mass sphere matrix with zeros so that each mass sphere for each atom has the same size.

        Args:
            mass_spheres_for_each_atom (list): mass spheres for each atom in molecule
        Returns:
            list: _description_
        """

        max_number_of_spheres = 0

        for spheres in mass_spheres_for_each_atom:

            max_number_of_spheres = max(max_number_of_spheres, len(spheres))

        for spheres in mass_spheres_for_each_atom:

            spheres += [[0]] * (max_number_of_spheres - len(spheres))

        for n in range(max_number_of_spheres):

            max_number_of_sphere_atoms = 0

            for j in range(self.n_atoms):

                max_number_of_sphere_atoms = max(
                    max_number_of_sphere_atoms, len(mass_spheres_for_each_atom[j][n])
                )

            for j in range(self.n_atoms):

                new_sphere = mass_spheres_for_each_atom[j][n] + [0] * (
                    max_number_of_sphere_atoms - len(mass_spheres_for_each_atom[j][n])
                )

                mass_spheres_for_each_atom[j][n] = new_sphere

        final_mass_spheres = []

        for j in range(self.n_atoms):

            final_mass_spheres.append(
                [atom for sphere in mass_spheres_for_each_atom[j] for atom in sphere]
            )

        return final_mass_spheres
