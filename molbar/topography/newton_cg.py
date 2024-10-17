import numpy as np
import os
from scipy.optimize import minimize
import molbar.fortranlib.analytical_derivatives as fder
import molbar.fortranlib.energy as fen
from molbar.helper.printer import Printer


class NetwonCG:

    def __init__(self, write_trj=False) -> None:
        """
        Initializes instance variables for storing idealization results
        """
        if not hasattr(self, "write_trj"):

            self.write_trj = write_trj

        # If true, the idealization trajectory will be written to a file
        if self.write_trj:

            self.opt_path = os.path.join(self.debug_path, f"{self.optfilename}.log")

            if os.path.isfile(self.opt_path):

                os.remove(self.opt_path)

    def run_idealization(self):

        # Flatten geometry array for optimization (needs to be 1D array)
        x0 = np.asfortranarray(self.coordinates.flatten())

        # Arguments for optimization such as bond, angle, dihedral parameters etc.
        self.args = (
            self.n_atoms,
            self.bonds,
            self.ideal_bonds,
            self.angles,
            self.ideal_angles,
            self.dihedrals,
            self.ideal_dihedrals,
            self.repulsion,
            self.charges,
            self.k_bond,
            self.k_angle,
            self.k_dihedral,
        )

        # Set options for optimization
        options = {"maxiter": 1e3, "xtol": 1e-15, "disp": False}

        # Run optimization
        # If write_trj is true, the trajectory will be written to a file
        if self.write_trj and self.debug_path is not None:

            res = minimize(
                self.get_energy,
                x0,
                method="Newton-CG",
                jac=self.get_gradient,
                hess=self.get_hessian,
                args=self.args,
                options=options,
                callback=self.print_opt_path,
            )

        else:

            res = minimize(
                self.get_energy,
                x0,
                method="Newton-CG",
                jac=self.get_gradient,
                hess=self.get_hessian,
                args=self.args,
                options=options,
            )

        # Store opt geometry and energy
        self.optimized_geometry = res.x.reshape(-1, 3)

        self.energy = res.fun

        self.convergence_type = res.success

    def print_opt_path(self, log):
        """
        Prints the idealization trajectory and energy to a file
        """
        # Reshape log array to 2D array (n_atoms, 3)
        current_geometry = log.reshape(-1, 3)

        # Get energy for geometry of trajectory
        energy = self.get_energy(
            log,
            self.n_atoms,
            self.bonds,
            self.ideal_bonds,
            self.angles,
            self.ideal_angles,
            self.dihedrals,
            self.ideal_dihedrals,
            self.repulsion,
            self.charges,
            self.k_bond,
            self.k_angle,
            self.k_dihedral,
        )

        # Print trajectory and energy to file

        Printer(
            self.n_atoms,
            energy,
            current_geometry,
            self.atoms,
            self.opt_path,
            update=True,
        ).print()

    def get_energy(
        self,
        x: np.ndarray,
        n_atoms: int,
        bonds: np.ndarray,
        ideal_bonds: np.ndarray,
        angles: np.ndarray,
        ideal_angles: np.ndarray,
        dihedrals: np.ndarray,
        ideal_dihedrals: np.ndarray,
        repulsion: np.ndarray,
        charges: np.ndarray,
        k_bond: float,
        k_angle: float,
        k_dihedral: float,
    ):
        """
        Calculates the energy of the current geometry

        Args:
            x (_type_): flattened geometry array
            n_atoms (int): number of atoms
            bonds (np.ndarray): atom indices of bonds
            ideal_bonds (np.ndarray): ideal bond lengths
            angles (np.ndarray): atom indices of anglesq
            ideal_angles (np.ndarray): ideal angle values
            dihedrals (np.ndarray): atom indices of dihedrals
            ideal_dihedrals (np.ndarray): ideal dihedral values
            repulsion (np.ndarray): atom indices of repulsion
            charges (np.ndarray): charges of atoms
            k_bond (float): force constant for bond
            k_angle (float): force constant for angle
            k_dihedral (float): force constant for dihedral

        Returns:
            energy (float): energy of current geometry
        """
        current_geometry = np.asfortranarray(x.reshape(-1, 3)).transpose()

        energy = fen.force_field_energy.get_energy(
            n_atoms,
            current_geometry,
            bonds,
            ideal_bonds,
            angles,
            ideal_angles,
            dihedrals,
            ideal_dihedrals,
            repulsion,
            charges,
            k_bond,
            k_angle,
            k_dihedral,
        )

        return energy

    def get_gradient(
        self,
        x: np.ndarray,
        n_atoms: int,
        bonds: np.ndarray,
        ideal_bonds: np.ndarray,
        angles: np.ndarray,
        ideal_angles: np.ndarray,
        dihedrals: np.ndarray,
        ideal_dihedrals: np.ndarray,
        repulsion: np.ndarray,
        charges: np.ndarray,
        k_bond: float,
        k_angle: float,
        k_dihedral: float,
    ):
        """
        Calculates the gradient of the current geometry

        Args:
            x (_type_): flattened geometry array
            n_atoms (int): number of atoms
            bonds (np.ndarray): atom indices of bonds
            ideal_bonds (np.ndarray): ideal bond lengths
            angles (np.ndarray): atom indices of anglesq
            ideal_angles (np.ndarray): ideal angle values
            dihedrals (np.ndarray): atom indices of dihedrals
            ideal_dihedrals (np.ndarray): ideal dihedral values
            repulsion (np.ndarray): atom indices of repulsion
            charges (np.ndarray): charges of atoms
            k_bond (float): force constant for bond
            k_angle (float): force constant for angle
            k_dihedral (float): force constant for dihedral

        Returns:
            energy (float): energy of current geometry
        """

        current_geometry = np.asfortranarray(x.reshape(-1, 3)).transpose()

        gradient, self.hessian = fder.derivatives.get_derivatives(
            n_atoms,
            current_geometry,
            bonds,
            ideal_bonds,
            angles,
            ideal_angles,
            dihedrals,
            ideal_dihedrals,
            repulsion,
            charges,
            k_bond,
            k_angle,
            k_dihedral,
        )

        return gradient

    def get_hessian(
        self,
        x: np.ndarray,
        n_atoms: int,
        bonds: np.ndarray,
        ideal_bonds: np.ndarray,
        angles: np.ndarray,
        ideal_angles: np.ndarray,
        dihedrals: np.ndarray,
        ideal_dihedrals: np.ndarray,
        repulsion: np.ndarray,
        charges: np.ndarray,
        k_bond: float,
        k_angle: float,
        k_dihedral: float,
    ):
        """
        Calculates the hessian of the current geometry

        Args:
            x (np.nadarray): flattened geometry array
            n_atoms (int): number of atoms
            bonds (np.ndarray): atom indices of bonds
            ideal_bonds (np.ndarray): ideal bond lengths
            angles (np.ndarray): atom indices of anglesq
            ideal_angles (np.ndarray): ideal angle values
            dihedrals (np.ndarray): atom indices of dihedrals
            ideal_dihedrals (np.ndarray): ideal dihedral values
            repulsion (np.ndarray): atom indices of repulsion
            charges (np.ndarray): charges of atoms
            k_bond (float): force constant for bond
            k_angle (float): force constant for angle
            k_dihedral (float): force constant for dihedral

        Returns:
            hessian (float): energy of current geometry
        """

        # current_geometry = np.asfortranarray(x.reshape(-1,3)).transpose()

        # hessian = fder.derivatives.get_hessian(n_atoms, current_geometry, bonds, ideal_bonds, angles, ideal_angles, dihedrals, ideal_dihedrals, repulsion, charges, k_bond, k_angle, k_dihedral)

        return self.hessian
