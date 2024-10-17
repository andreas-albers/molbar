import numpy as np


class Vector:
    """
    A class to represent a set of vectors.

    Args:
        _i (np.ndarray): first vector
        _j (np.ndarray): second vector
        _k (np.ndarray, optional): third vector. Defaults to None.
        _l (np.ndarray, optional): fourth vector. Defaults to None.

    Methods:

        unit_vector(): Calculates unit vector between two vectors.
        calculate_distance(): Calculates distance between two vectors i and j.
        calculate_angle(): Calculates the angle between three vectors i, j and k.
        calculate_dihedral(): Calculates the dihedral between four vectors i, j, k and l.
    """

    def __init__(self, cartesian_points) -> None:
        """
        Constructs all the necassary attributes for the Vector object.

        Args:
            cartesian_points (np.ndarray): set of cartesian points

        """

        self.cartesian_points = cartesian_points

    def distance_vector(self) -> np.ndarray:
        """
        Calculates distance vector between two vectors.

        Returns:
            distance vector (np.ndarray):  distance vector between the i and j vector

        """

        return self.cartesian_points[0] - self.cartesian_points[1]

    def unit_vector(self) -> np.ndarray:
        """
        Calculates unit vector between two vectors.

        Returns:
            unit vector (np.ndarray):  unit vector between two vectors.

        """

        distance_vector = self.distance_vector()

        if np.linalg.norm(distance_vector) != 0:

            return distance_vector / np.linalg.norm(distance_vector)

        else:

            return np.zeros(3)

    def calculate_distance(self) -> float:
        """
        Calculates distance between two vectors i and j.

        Returns:
            distance (float): distance between the two vectors.

        """

        distance_vector = self.distance_vector()

        return np.linalg.norm(distance_vector)

    def calculate_angle(self) -> float:
        """
        Calculates the angle between three vectors.

        Returns:
            distance (float): angle between the i,j and k vector
        """

        r_ij = self.cartesian_points[0] - self.cartesian_points[1]

        r_kj = self.cartesian_points[2] - self.cartesian_points[1]

        cos_theta = np.dot(r_ij, r_kj)

        sin_theta = np.linalg.norm(np.cross(r_ij, r_kj))

        theta = np.arctan2(sin_theta, cos_theta)

        return 180.0 * theta / np.pi

    def calculate_dihedral(self) -> float:
        """
        Calculates the dihedral between four vectors.

        Returns:
            float: dihedral between i, j, k and l vector
        """

        b = self.cartesian_points[:-1] - self.cartesian_points[1:]

        b[0] *= -1

        v = np.array([v - (v.dot(b[1]) / b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])

        # Normalize vectors

        v /= np.sqrt(np.einsum("...i,...i", v, v)).reshape(-1, 1)

        b1 = b[1] / np.linalg.norm(b[1])

        x = np.dot(v[0], v[1])

        m = np.cross(v[0], b1)

        y = np.dot(m, v[1])

        return np.degrees(np.arctan2(y, x))
