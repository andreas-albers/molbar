import numpy as np


class Printer:

    def __init__(
        self,
        n_atoms: int,
        energy: float,
        coordinates: np.ndarray,
        elements: np.ndarray,
        path: str,
        update=False,
    ) -> None:
        """
        Constructs a Printer object.

        Args:
            n_atoms (int): number of atoms
            energy (float): energy of the structure
            coordinates (np.ndarray): cartesian coordinates of the structure
            elements (np.ndarray): elements of the structure
            path (str): path to the output file
            update (bool, optional): Whether to append existing .xyz file. Defaults to False.
        """

        self.n_atoms = n_atoms

        self.energy = energy

        self.coordinates = coordinates

        self.elements = elements

        self.path = path

        if update == True:

            self.edit_mode = "a"

        else:

            self.edit_mode = "w"

    def print(self):

        with open(self.path, self.edit_mode) as output:

            print(self.n_atoms, file=output)
            print(f"energy: {self.energy}", file=output)

            for j in range(self.n_atoms):

                print(
                    "{:<4s}\t{:>11.12f}\t{:>11.12f}\t{:>11.12f}".format(
                        self.elements[j],
                        self.coordinates[j][0],
                        self.coordinates[j][1],
                        self.coordinates[j][2],
                    ),
                    file=output,
                )
