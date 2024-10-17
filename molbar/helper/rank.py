import numpy as np


class Rank:
    """
    A class for representing the prioritization of substituents around a central atom.

    Args:

        _rank_matrix (list): rank matrix, ranking the nth mass sphere of the ith group with other nth mass sphere with other j groups.

    Method:
        get_ranking(): Identifiying uniques in rank list.

    """

    def __init__(self, rank_matrix: np.array) -> None:
        """
        Constructs all the necassary attributes for the Rank object.
        """

        self.rank_matrix = rank_matrix

        self.get_ranking()

    def get_ranking(self) -> np.ndarray:
        """
        Map 2D-rank matrix to final 1D-rank for each substituent.

        Returns:
            np.ndarray: final 1D-rank for each substituent.
        """

        number_of_spheres = len(self.rank_matrix[0])

        keys = [str(i) for i in range(number_of_spheres)]

        dtype = [(str(i), int) for i in range(number_of_spheres)]

        values = [tuple(molecule) for molecule in self.rank_matrix]

        rank = np.array(values, dtype)

        sorted_index = np.argsort(rank, order=keys, axis=0, kind="mergesort")

        rank_without_duplicates = [
            np.where(sorted_index == molecule)[0][0]
            for molecule in range(len(sorted_index))
        ]

        _, starting_index, counts = np.unique(
            rank, return_counts=True, return_index=True
        )

        final_rank = rank_without_duplicates.copy()

        for index, count in zip(starting_index, counts):

            rank_for_duplicates = rank_without_duplicates[index]

            for i in range(1, count):

                duplicate_index = np.where(
                    rank_without_duplicates == (rank_for_duplicates + i)
                )[0][0]

                final_rank[duplicate_index] = rank_for_duplicates

        self.priorities = np.array(final_rank) + 1
