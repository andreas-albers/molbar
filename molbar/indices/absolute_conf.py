import numpy as np
import molbar.fortranlib.opd as opd
from molbar.helper.symcheck import SymCheck
from molbar.helper.printer import Printer
from scipy.spatial.transform import Rotation


class AbsoluteConfIndex:

    def __init__(self, optimized_geometry, priorities, core_indices):
        """
        Initializes instance variables for idealization.

        Parameters
        ----------
        optimized_geometry : np.ndarray
            The optimized geometry of the molecule.
        priorities : np.ndarray
            The priorities of the atoms.
        core_indices : np.ndarray
            The indices of the core atoms.

        Returns
        -------
        None
        """

        self.optimized_geometry = optimized_geometry

        self.priorities = priorities

        self.core_indices = core_indices

        self.n_atoms = len(optimized_geometry)

        self.n_core = len(core_indices)

    def get_absconf_index(self) -> int:
        """
        Determines absolute configuration index based on "M. Osipov, B. Pickup, D. Dunmur, Mol. Phys. 1995, 84, 1193–1206."

        Returns:
            int: chirality index

        """
        if self.n_atoms == 4 and self.n_core == 1:
            return self.get_bannwarth_3()
        elif self.n_atoms == 5 and self.n_core == 1:
            return self.get_bannwarth_4()
        elif self.n_atoms > 2:
            return self.get_odp()
        else:
            return 0

    def get_odp(self) -> int:
        """
        Gets absolute configuration index based on the ODP index  "M. Osipov, B. Pickup, D. Dunmur, Mol. Phys. 1995, 84, 1193–1206."

        Returns:
            int: ODP index
        """

        # Transpose the final idealized geometry array to be in Fortran order
        fgeometry = np.asfortranarray(self.optimized_geometry).transpose()

        # Convert the priorities array to Fortran order
        # fpriorities = np.asfortranarray(rankdata(self.priorities, method="min"))

        fpriorities = np.asfortranarray(self.priorities)
        if not self.check_chiral():
            return 0
        # Get the ODP index from fortran extension
        G_os = opd.opd.get_index(self.n_atoms, fgeometry, fpriorities)
        # G_os = sel.get_opd(self.n_atoms, sorted_coordinates, sorted_priorities)
        if G_os > 1e-6:
            return 1
        elif G_os < -1e-6:
            return -1
        else:
            return 0

    def get_opd_python(self, n_atoms, coordinates, priorities):

        G = 0.0

        for i in range(n_atoms):
            for j in range(n_atoms):
                for k in range(n_atoms):
                    for l in range(n_atoms):

                        if i != j and k != l and i != l and j != k:

                            G += self.get_gijkl(i, j, k, l, coordinates, priorities)

        return 24.0 / (float(n_atoms) ** 4) * G
        # return G

    def get_gijkl(self, i, j, k, l, coordinates, priorities):

        pi = priorities[i]
        pj = priorities[j]
        pk = priorities[k]
        pl = priorities[l]

        ri = coordinates[i]
        rj = coordinates[j]
        rk = coordinates[k]
        rl = coordinates[l]
        rij = ri - rj
        rkl = rk - rl
        ril = ri - rl
        rjk = rj - rk

        norm_rij = np.linalg.norm(rij)
        norm_rjk = np.linalg.norm(rjk)
        norm_rkl = np.linalg.norm(rkl)
        norm_ril = np.linalg.norm(ril)

        # print(rij, norm_rij)
        # print(rkl, norm_rkl)
        # print(rjk, norm_rjk)
        # print(ril, norm_ril)

        cross_ijkl = np.cross(rij, rkl)

        dot_ijjk = np.dot(rij, rjk)

        dot_jkkl = np.dot(rjk, rkl)

        # print(cross_ijkl)
        # print(dot_ijjk)
        # print(dot_jkkl)
        # print(np.dot(cross_ijkl, ril))
        # print(pi*pj*pk*pl)

        enumerator = np.dot(cross_ijkl, ril) * dot_ijjk * dot_jkkl

        denominator = (norm_rij * norm_rjk * norm_rkl) ** 2 * norm_ril

        return pi * pj * pk * pl * enumerator / denominator

    def get_bannwarth_4(self):
        """
        Determines the chirality index for a tetrahedral structure based on the triple product of the three highest prioritized adjacent atoms.

        Returns:
            int: chirality index
        """

        core_index = self.core_indices[0]

        # Define the geometry array as the difference between the idealized geometry and the core atom
        geometry = np.array(
            [
                self.optimized_geometry[i] - self.optimized_geometry[core_index]
                for i in range(self.n_atoms)
                if i != core_index
            ]
        )

        # Define the priorities array as the priorities of the adjacent atoms
        prio = np.array(
            [self.priorities[i] for i in range(self.n_atoms) if i != core_index]
        )

        # If there are four unique priorities, calculate the determinant of the geometry array
        if len(set(prio)) == 4:
            # Get the indices of the three highest priorities
            idx = prio.argsort()[-3:]
            # Calculate the determinant of the geometry array
            det = np.linalg.det(
                np.array([geometry[idx[0]], geometry[idx[1]], geometry[idx[2]]])
            )

            if det > 0.1:

                return 1

            elif det < 0.1:

                return -1

            else:

                return 0

        else:

            return 0

    def get_bannwarth_3(self):
        """
        Determines the chirality index for a trigonal pyramidal structure based on the triple product of the three highest prioritized adjacent atoms.

        Returns:
            int: chirality index
        """
        core_index = self.core_indices[0]

        # Define the geometry array as the difference between the idealized geometry and the core atom
        geometry = np.array(
            [
                self.optimized_geometry[i] - self.optimized_geometry[core_index]
                for i in range(self.n_atoms)
                if i != core_index
            ]
        )

        # Define the priorities array as the priorities of the adjacent atoms
        prio = np.array(
            [self.priorities[i] for i in range(self.n_atoms) if i != core_index]
        )
        # If there are three unique priorities, calculate the determinant of the geometry array
        if len(set(prio)) == 3:
            # Get the indices of the three highest priorities
            idx = prio.argsort()[-3:]
            # Calculate the determinant of the geometry array
            det = np.linalg.det(
                np.array([geometry[idx[0]], geometry[idx[1]], geometry[idx[2]]])
            )

            if det > 0.5:

                return 1

            elif det < -0.5:

                return -1

            else:

                return 0

        else:

            return 0

    def check_chiral(self):
        symCheck = SymCheck(self.optimized_geometry, self.priorities)
        if not symCheck.check_for_possible_chiral():
            return False
        centered_coordinates = symCheck.return_aligned_coords()
        centered_coordinates_inverted = -1 * centered_coordinates.copy()
        R, _ = Rotation.align_vectors(
            centered_coordinates, centered_coordinates_inverted
        )
        coordinates_inverted_rotated = R.apply(centered_coordinates_inverted)
        rmsd = np.sqrt(
            np.sum((centered_coordinates - coordinates_inverted_rotated) ** 2)
            / len(centered_coordinates)
        )
        if rmsd < 0.3:
            return False
        else:
            return True
