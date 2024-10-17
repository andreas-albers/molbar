import numpy as np
from molbar.topography.newton_cg import NetwonCG


class GeomOpt(NetwonCG):

    def __init__(
        self,
        n_atoms=None,
        n_core=None,
        core_indices=None,
        coordinates=None,
        atoms=None,
        bonds=None,
        ideal_bonds=None,
        angles=None,
        ideal_angles=None,
        dihedrals=None,
        ideal_dihedrals=None,
        repulsion=None,
        charges=None,
        k_bond=1e6,
        k_angle=2e3,
        k_dihedral=2e3,
        write_trj=False,
        get_index=False,
        fragment_id=0,
        debug_path="debug",
    ) -> None:
        """
        Initializes instance variables for idealization.

        Parameters
        ----------
        write_trj : bool, optional
            Whether to write trajectory file, by default False
        get_index : bool, optional
            Whether to get relative and absolute configuration indices, by default False

        Returns
        -------
        None
        """

        # Set basic parameters

        self.n_atoms = n_atoms

        self.n_core = n_core

        self.core_indices = core_indices

        self.coordinates = coordinates

        self.atoms = atoms

        # Set constraints
        if (bonds is not None) and (len(bonds) > 0):

            self.bonds = np.asfortranarray(
                [[bond[0] + 1, bond[1] + 1] for bond in bonds]
            ).T

        else:

            self.bonds = np.asfortranarray([])

        self.ideal_bonds = np.asfortranarray(ideal_bonds)

        if (angles is not None) and (len(angles) > 0):

            self.angles = np.asfortranarray(
                [[angle[0] + 1, angle[1] + 1, angle[2] + 1] for angle in angles]
            ).T
        else:

            self.angles = np.asfortranarray([])

        if (ideal_angles is not None) and (len(ideal_angles) > 0):

            self.ideal_angles = np.asfortranarray(ideal_angles)

        else:

            self.ideal_angles = np.asfortranarray([])

        if (dihedrals is not None) and (len(dihedrals) > 0):

            self.dihedrals = np.asfortranarray(
                [
                    [dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1]
                    for dihedral in dihedrals
                ]
            ).T

        else:

            self.dihedrals = np.asfortranarray([])

        if (ideal_dihedrals is not None) and (len(ideal_dihedrals) > 0):

            self.ideal_dihedrals = np.asfortranarray(ideal_dihedrals)

        else:

            self.ideal_dihedrals = np.asfortranarray([])

        if (repulsion is not None) and len(repulsion) > 0:

            self.repulsion = np.asfortranarray(
                [[repulsive[0] + 1, repulsive[1] + 1] for repulsive in repulsion]
            ).T

        else:

            self.repulsion = np.asfortranarray([])

        if (charges is not None) and len(charges) > 0:

            self.charges = np.asfortranarray(charges)

        else:

            self.charges = np.asfortranarray([])

        self.fragment_id = fragment_id

        self.debug_path = debug_path

        # Initialize instance variables for storing idealization results

        self.write_trj = write_trj

        self.get_index = get_index

        # Force constants for bond, angle, and dihedral constraints

        self.k_bond = k_bond

        self.k_angle = k_angle

        self.k_dihedral = k_dihedral

    def optimize(self):
        """
        Idealizes geometry with force field.

        Returns
        -------
        None
        """

        self.check_parameters()

        # Set Idealization object parameters

        self.optfilename = str(self.fragment_id + 1)

        # Initialize NetwonCG object
        NetwonCG.__init__(self)

        # Run Newton-CG optimization with force field
        self.run_idealization()

    def check_parameters(self):
        """
        Checks if all necessary parameters are provided for force field initialization.

        Raises:
            ValueError: If any of the necessary parameters are missing.
        """

        necessary_parameters = np.array(
            [
                "n_atoms",
                "coordinates",
                "atoms",
                "bonds",
                "ideal_bonds",
                "angles",
                "ideal_angles",
                "dihedrals",
                "ideal_dihedrals",
                "repulsion",
                "charges",
            ]
        )

        parameters = vars(self)

        checked = np.array([parameters[key] is None for key in necessary_parameters])

        if any(checked):

            raise ValueError(
                "The following parameters are missing for force field initalization: {}".format(
                    necessary_parameters[np.where(checked == True)]
                )
            )

        if len(self.bonds) > 0:
            if self.bonds.shape[1] != len(self.ideal_bonds):

                raise ValueError(
                    "Bond and ideal bond arrays have different number of columns"
                )

        if len(self.angles) > 0:

            if self.angles.shape[1] != len(self.ideal_angles):

                raise ValueError(
                    "Angle and ideal angle arrays have different number of columns"
                )

        if len(self.dihedrals) > 0:
            if self.dihedrals.shape[1] != len(self.ideal_dihedrals):

                raise ValueError(
                    "Dihedral and ideal dihedral arrays have different number of columns"
                )

        if len(self.repulsion) > 0:
            if self.repulsion.shape[1] != len(self.charges):

                raise ValueError(
                    "Repulsion and charge arrays have different number of columns"
                )
