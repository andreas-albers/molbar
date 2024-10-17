import numpy as np
from ase.io import read


class FileReader:
    """
    A class to read a file containing molecular information such as coordinates and elements.

    Methods:
        read_file(): Reads number of atoms in the molecule and geometry from either a .xyz-, coord- or .sdf/.mol-file.

    """

    def __init__(self, filepath=None):
        """
        Initializes the FileReader class.
        Args:
            filepath (str): Path to file.
        """

        if not hasattr(self, "filepath"):

            self.filepath = filepath

    def read_file(self, filepath="") -> tuple:
        """
        Reads number of atoms in the molecule, cartesian coordinates and present atoms from .xyz file.

        Args:
            filepath (str): Path to file. Defaults to self.filepath.
            The file can be in formats such as XYZ, Turbomole coord, SDF/MOL, CIF, PDB, Gaussian (.com, .gjf, .log), and many others. For a complete list, refer to the ASE documentation, ASE I/O Formats: https://wiki.fysik.dtu.dk/ase/ase/io/io.html).
        Returns:
            Tuple containing the number of atoms in the molecule (int), cartesian coordinates (np.ndarray) and present atoms (np.ndarray).

        Raises:
            XYZNotFound: If the file does not exist.
            NotXYZFormat: If the input geometry is not in .xyz format.
        """

        molecule = read(self.filepath)
        # Get the atomic symbols and positions
        atoms = (
            molecule.get_chemical_symbols()
        )  # This gives you a list of atomic symbols
        geometry = np.array(molecule.get_positions())
        n_atoms = len(atoms)
        return n_atoms, geometry, atoms
