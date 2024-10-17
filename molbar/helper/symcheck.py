import numpy as np
from scipy.optimize import linear_sum_assignment, minimize
from scipy.spatial.transform import Rotation as R


class SymCheck:
    """
    A class to check for symmetry operations in a molecule.

    Attributes:
    coordinates (np.ndarray): The coordinates of the molecule.
    n_atoms (int): The number of atoms in the molecule.
    masses (np.ndarray): The masses of the atoms.
    threshold (float): The threshold for the cost function.
    rotated_coords (np.ndarray): The aligned coordinates.
    moments (np.ndarray): The moments of inertia.
    principal_axes (np.ndarray): The principal axes.
    mask (np.ndarray): The mask to prevent the same atom from being matched non-equivalent atoms in the Hungarian algorithm.
    """

    def __init__(self, coordinates, masses, threshold=1e-3) -> None:
        self.coordinates = coordinates
        self.n_atoms = len(coordinates)
        self.masses = masses
        self.threshold = threshold
        self.rotated_coords, self.moments, self.principal_axes = (
            self.align_to_principal_axes()
        )
        self.mask = self.get_mask()

    def return_aligned_coords(self):
        """
        Align the coordinates to the principal axes.
        Returns:
        np.ndarray: The aligned coordinates.
        """
        return self.rotated_coords

    def return_moments_of_inertia(self):
        """
        Return the moments of inertia.
        Returns:
        np.ndarray: The moments of inertia.
        """
        return self.moments

    def return_principal_axes(self):
        """
        Return the principal axes.
        Returns:
        np.ndarray: The principal axes.
        """
        return self.principal_axes

    def return_center_of_mass(self):
        """
        Calculate the center of mass of the molecule.
        Returns:
        np.ndarray: The center of mass.
        """
        total_mass = np.sum(self.masses)
        return np.sum(self.coordinates.T * self.masses, axis=1) / total_mass

    def translate_to_origin(self):
        """
        Translate the molecule to the origin.
        Returns:
        np.ndarray: The translated coordinates.
        """
        com = self.return_center_of_mass()
        return self.coordinates - com

    def get_inertia_tensor(self, coordinates, masses):
        """
        Calculate the inertia tensor of the molecule.
        Parameters:
        coordinates (np.ndarray): The translated coordinates.
        masses (np.ndarray): The masses of the atoms.
        Returns:
        np.ndarray: The inertia tensor.
        """
        Ixx = np.sum(masses * (coordinates[:, 1] ** 2 + coordinates[:, 2] ** 2))
        Iyy = np.sum(masses * (coordinates[:, 0] ** 2 + coordinates[:, 2] ** 2))
        Izz = np.sum(masses * (coordinates[:, 0] ** 2 + coordinates[:, 1] ** 2))
        Ixy = -np.sum(masses * coordinates[:, 0] * coordinates[:, 1])
        Ixz = -np.sum(masses * coordinates[:, 0] * coordinates[:, 2])
        Iyz = -np.sum(masses * coordinates[:, 1] * coordinates[:, 2])
        inertia_tensor = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
        return inertia_tensor

    def align_to_principal_axes(self):
        """
        Align the molecule to the principal axes.
        Returns:
        np.ndarray: The aligned coordinates.
        """
        translated_coords = self.translate_to_origin()
        I_tensor = self.get_inertia_tensor(translated_coords, self.masses)
        eigenvalues, eigenvectors = np.linalg.eig(I_tensor)
        idx = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        rotated_coords = np.dot(translated_coords, eigenvectors)
        return rotated_coords, eigenvalues, eigenvectors

    def get_mask(self):
        """
        Create a mask to prevent the same atom from being matched non-equivalent atoms in the Hungarian algorithm when checking if two atoms are on the same positions before and after a symmetry operation.
        Returns:
        np.ndarray: The mask.
        """
        unique_masses = np.unique(self.masses)
        mask = np.ones((self.n_atoms, self.n_atoms, len(unique_masses)), dtype=bool)
        for idx, mass in enumerate(unique_masses):
            group_indices = np.where(self.masses == mass)[0]
            mask_layer = np.ones((self.n_atoms, self.n_atoms), dtype=bool)
            mask_layer[np.ix_(group_indices, group_indices)] = False
            mask[:, :, idx] = mask_layer
        final_mask = np.all(mask, axis=-1)
        return final_mask

    def get_cost_function(self, coordinates, flipped_coordinates):
        """
        Calculate the cost function for the Hungarian algorithm.
        Parameters:
        coordinates (np.ndarray): The original coordinates.
        flipped_coordinates (np.ndarray): The flipped coordinates.
        Returns:
        np.ndarray: The cost matrix.
        """
        large_value = np.inf
        cost_matrix = np.abs(flipped_coordinates[:, np.newaxis] - coordinates).sum(
            axis=2
        )
        cost_matrix[self.mask] = large_value
        return cost_matrix

    def get_cost_for_symmetry_operation(self, operation):
        """
        Check if a symmetry operation is present in the molecule.
        Parameters:
        operation (str): The symmetry operation to check for. Either "xy", "yz", "xz", or "inv".
        threshold (float): The threshold for the cost function.
        Returns:
        bool: True if the symmetry operation is present, False otherwise.
        """
        flipped_coordinates = np.copy(self.rotated_coords)
        if operation == "xy":
            flipped_coordinates[:, 2] *= -1
        elif operation == "yz":
            flipped_coordinates[:, 0] *= -1
        elif operation == "xz":
            flipped_coordinates[:, 1] *= -1
        elif operation == "inv":
            flipped_coordinates *= -1
        else:
            raise ValueError(f"Invalid operation: {operation}")
        cost_matrix = self.get_cost_function(self.rotated_coords, flipped_coordinates)
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        total_cost = cost_matrix[row_ind, col_ind].sum().mean()
        return total_cost

    def check_for_possible_chiral(self):
        """
        Check if the molecule is chiral.
        Parameters:
        threshold (float): The threshold for the cost function.
        Returns:
        bool: True if the molecule is possible chiral, False otherwise.
        """
        n_distinct = self.check_for_distinct_moments()
        if n_distinct == 0:
            return False
        elif n_distinct == 1:
            if self.handle_non_asymmetric() < self.threshold:
                return False
            return True
        operations = ["xy", "yz", "xz", "inv"]
        for operation in operations:
            if self.get_cost_for_symmetry_operation(operation) < self.threshold:
                return False
        return True

    def handle_non_asymmetric(self):
        if self.distinct_moments[0] == 0:
            axis = "x"
            symmetry_operations = ["xy", "xz"]
        elif self.distinct_moments[0] == 1:
            axis = "y"
            symmetry_operations = ["xy", "yz"]
        else:
            axis = "z"
            symmetry_operations = ["xz", "yz"]
        costs = []
        for operation in symmetry_operations:
            cost = self.get_cost_for_symmetry_operation(operation)
            costs.append(cost)
        best_symmetry_operation = symmetry_operations[np.argmin(costs)]
        self.rotated_coords_backup = np.copy(self.rotated_coords)

        def objective_function(angle, axis, best_symmetry_operation):
            angle_radians = np.deg2rad(angle)
            # Create a rotation object for rotation around the z-axis
            rotation = R.from_euler(axis, angle_radians)
            # Apply the rotation to the coordinates
            self.rotated_coords = rotation.apply(self.rotated_coords_backup)
            cost = self.get_cost_for_symmetry_operation(best_symmetry_operation)
            return cost

        result = minimize(
            objective_function,
            0.0,
            method="Nelder-Mead",
            args=(axis, best_symmetry_operation),
        )
        return result.fun

    def check_for_distinct_moments(self, threshold=1e-1):
        """
        Check if all values in a NumPy array are distinct within a certain threshold.
        Parameters:
        arr (np.ndarray): Input array to check.
        threshold (float): The minimum difference between any two values to consider them distinct.
        Returns:
        bool: True if all values are distinct within the threshold, False otherwise.
        """
        differences = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                if j > i:
                    moments = [self.moments[i], self.moments[j]]
                    if np.abs(np.diff(moments)) < threshold:
                        differences[i, j] = 1
                        differences[j, i] = 1
        cost = np.sum(differences, axis=0)
        n_distinct = np.sum(cost == 0)
        self.distinct_moments = np.where(cost == 0)[0]
        return n_distinct


def main():
    import sys
    from molbar.molgraph.molgraph import MolGraph
    from molbar.helper.printer import Printer
    import os
    import argparse

    # Step 1: Create the parser
    parser = argparse.ArgumentParser(
        description="Helper script to split an align the molecule to the principal axes."
    )

    # Step 2: Add arguments
    parser.add_argument("file", type=str, help="Input file")
    parser.add_argument(
        "-r",
        "--replace",
        action="store_true",
        help="Overwrite the input file with the aligned coordinates.",
    )
    # Step 3: Parse the arguments
    args = parser.parse_args()
    path = args.file
    if args.replace:
        new_path = path
    else:
        directory = os.path.dirname(path)
        filename = os.path.basename(path)
        new_filename = filename.split(".")[0] + "_aligned.xyz"
        new_path = os.path.join(directory, new_filename)
    molGraph = MolGraph()
    molGraph.from_file(path)
    coordinates = molGraph.return_node_data(attribute="coordinates")
    masses = molGraph.return_node_data(attribute="masses")
    symcheck = SymCheck(coordinates, masses)
    final_geometry = symcheck.return_aligned_coords()
    Printer(
        n_atoms=molGraph.return_n_atoms(),
        energy="",
        coordinates=final_geometry,
        elements=molGraph.return_node_data(attribute="elements"),
        path=new_path,
    ).print()


if __name__ == "__main__":
    from molbar.io.filereader import FileReader
    from molbar.helper.printer import Printer
    import json

    name = "c76_m"
    path = f""
    idx = 5
    reader = FileReader(f"{path}{idx}.xyz")
    n_atoms, coordinates, elements = reader.read_file()
    with open(f"{path}/output.json", "r") as f:
        data = json.load(f)
        masses = list(data["fragment_data"][str(idx)]["priorities"].values())
    symcheck = SymCheck(coordinates, masses)
    final_geometry = symcheck.return_aligned_coords()
    possible_chiral = symcheck.check_for_possible_chiral()
    Printer(
        n_atoms=n_atoms,
        energy="",
        coordinates=final_geometry,
        elements=elements,
        path=f"{path}trafo.xyz",
    ).print()
