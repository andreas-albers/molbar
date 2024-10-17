import numpy as np
import os
from molbar.barcode import get_molbars_from_coordinates
from molbar.helper.printer import Printer


class EnsembleSplitter:

    def __init__(
        self, filename=None, list_of_coordinates=None, list_of_elements=None
    ) -> None:
        assert (
            filename is not None or list_of_coordinates is not None
        ), "Either a filename or a list of coordinates must be provided."
        if filename is not None:
            assert (
                list_of_coordinates is None
            ), "If a filename is provided, the list of coordinates must be None."
            assert isinstance(filename, str), "The filename must be a string."
            self.filename = filename
            self.list_of_coordinates, self.list_of_elements, self.list_of_comments = (
                self.split_file_ensemble(filename)
            )
        if list_of_coordinates is not None:
            assert (
                list_of_elements is not None
            ), "If a list of coordinates is provided, the list of elements must be provided."
            assert (
                filename is None
            ), "If a list of coordinates is provided, the filename must be None."
            assert isinstance(
                list_of_elements, np.ndarray
            ), "The list of elements must be a numpy array."
            assert len(list_of_coordinates) == len(
                list_of_elements
            ), "The number of molecules provided by elements or by coordinates must be the same."
            self.list_of_coordinates = list_of_coordinates
            self.list_of_elements = list_of_elements
            self.list_of_comments = None

    def split_ensemble(self, threads=1, show_progress_bar=False, topo=False):
        """
        Splits a collection of ensembles into single ensembles where each ensemble is characterized by a unique molecular barcode, so it containts only conformers.

        Args:
            threads (int, optional): The number of threads to use. Defaults to 1.
            show_progress_bar (bool, optional): Whether to show a progress bar. Defaults to False.
        """
        self.molbars = self.get_molbars(
            threads=threads, progress_bar=show_progress_bar, topo=topo
        )
        self.unique_molbars = np.unique(self.molbars)
        ensemble_idxs = self.order_ensembles()
        ensembles_with_molID = self.group_molecules_into_ensembles()
        ensembles_with_molID, boltzmann_averages = self.get_ensemble_energies(
            ensembles_with_molID
        )
        ordered_grouped_molecules = []
        file_names = []
        grouped_boltzman_averages = [
            np.mean([boltzmann_averages[idx] for idx in idx_group])
            for idx_group in ensemble_idxs
        ]
        argsort_ensemble_idxs = np.argsort(grouped_boltzman_averages)
        for ensemble_idx in argsort_ensemble_idxs:
            group = ensemble_idxs[ensemble_idx]
            file_name = f"ensemble_{ensemble_idx+1}"
            ordered_grouped_molecules.append(ensembles_with_molID[group[0]])
            file_names.append(file_name + ".xyz")
            if len(group) == 2:
                ordered_grouped_molecules.append(ensembles_with_molID[group[1]])
                file_names.append(file_name + "_enantiomer.xyz")
        self.ordered_grouped_molecules = ordered_grouped_molecules
        self.file_names = file_names

    def group_molecules_into_ensembles(self):
        """
        Groups the molecules into ensembles based on their molecular barcode.
        """
        groups = []
        for unique_molbar in self.unique_molbars:
            group = []
            for idx, molbar in enumerate(self.molbars):
                if molbar == unique_molbar:
                    group.append(idx)
            groups.append(group)
        return groups

    def get_ensemble_energies(self, groups):

        def boltzmann_average(energies, temperature):
            k_B = 1.380649e-23  # Boltzmann constant in J/K
            min_energy = np.min(energies)
            boltzmann_factors = np.exp(
                -(np.array(energies) - min_energy) / (k_B * temperature)
            )
            partition_function = np.sum(boltzmann_factors)
            average_energy = (
                np.sum(np.array(energies) * boltzmann_factors) / partition_function
            )
            return average_energy

        def is_float(value):
            try:
                float(value)
                return True
            except ValueError:
                return False

        sorted_groups = []
        boltzmann_averages = []
        for idxs in groups:
            ensemble_comments = [self.list_of_comments[idx] for idx in idxs]
            try:
                energies = []
                for comment in ensemble_comments:
                    split_comment = comment.split()
                    found_energy = False
                    for word in split_comment:
                        if is_float(word):
                            energies.append(float(word))
                            found_energy = True
                            break
                    if not found_energy:
                        energies.append(None)
                sorted_idxs = np.argsort(energies)
                new_idxs = np.array(idxs)[sorted_idxs]
                sorted_groups.append(new_idxs)
            except ValueError:
                energies = None
            boltzmann_averages.append(boltzmann_average(energies, 298.15))
        return sorted_groups, boltzmann_averages

    def return_ensembles(self):
        """
        Returns the ensembles as a list of lists where each list contains the coordinates, elements, and comments of a molecule.
        """
        list_of_ensembles = []
        for idxs in self.ordered_grouped_molecules:
            ensemble = []
            for idx in idxs:
                ensemble.append(
                    [
                        self.list_of_coordinates[idx],
                        self.list_of_elements[idx],
                        self.list_of_comments[idx],
                    ]
                )
            list_of_ensembles.append(ensemble)
        return list_of_ensembles

    def print_ensembles(self, target_directory: str):
        """
        Prints the ensembles to a directory.

        Args:
            target_directory (str): The target directory where the ensembles should be printed.
        """

        if not os.path.exists(target_directory):
            os.makedirs(target_directory)
        for idxs, filename in zip(self.ordered_grouped_molecules, self.file_names):
            full_path = os.path.join(target_directory, filename)
            if os.path.isfile(full_path):
                os.remove(full_path)

            ensemble_coordinates = [self.list_of_coordinates[idx] for idx in idxs]
            ensemble_elements = [self.list_of_elements[idx] for idx in idxs]
            ensemble_comments = [self.list_of_comments[idx] for idx in idxs]
            for idx in idxs:
                Printer(
                    n_atoms=len(self.list_of_coordinates[idx]),
                    energy=self.list_of_comments[idx],
                    coordinates=self.list_of_coordinates[idx],
                    elements=self.list_of_elements[idx],
                    path=full_path,
                    update=True,
                ).print()

    def order_ensembles(self):
        """
        Orders the ensembles based on the molecular barcode. Ensembles with the same topography and opposite chirality spectrum are considered enantiomers.

        Returns:
            list: A list of lists where each list contains the indices of the molecules in an ensemble.
        """
        splitted_molbars = [molbar.rsplit("|", 1) for molbar in self.unique_molbars]
        enantiomers = []
        single = []
        for i in range(len(self.unique_molbars)):
            if i in np.array(enantiomers).flatten():
                continue
            any_same_topography = False
            for j in range(i + 1, len(self.unique_molbars)):
                ith_topography = splitted_molbars[i][0]
                jth_topography = splitted_molbars[j][0]
                ith_chirality = np.array(
                    [int(x) for x in splitted_molbars[i][1].split()]
                )
                jth_chirality = sorted(
                    -1 * np.array([int(x) for x in splitted_molbars[j][1].split()])
                )
                if ith_topography == jth_topography:
                    if np.all(ith_chirality == jth_chirality):
                        enantiomers.append([i, j])
                        any_same_topography = True
            if not any_same_topography:
                single.append([i])
        order = enantiomers + single
        return order

    def get_molbars(self, threads=1, progress_bar=False, topo=False):
        """
        Returns the molecular barcodes for the ensemble.

        Args:
            threads (int, optional): The number of threads to use. Defaults to 1.
            progress_bar (bool, optional): Whether to show a progress bar. Defaults to False.
        """
        if topo:
            mode = "topo"
        else:
            mode = "mb"
        molbars = get_molbars_from_coordinates(
            self.list_of_coordinates,
            self.list_of_elements,
            threads=threads,
            progress=progress_bar,
            mode=mode,
        )
        return molbars

    def split_file_ensemble(self, filename):
        """
        Splits an ensemble file into coordinates, elements, and comments.

        Args:
            filename (str): The filename of the ensemble file.

        Returns:
            list: A list of coordinates.
            list: A list of elements.
            list: A list of comments.
        """
        with open(filename, "r") as f:
            lines = f.readlines()
            n_atoms = 0
            counter = -1
            list_of_coordinates = []
            list_of_elements = []
            list_of_comments = []
            coordinates = []
            elements = []
            counter = -1
            for idx, line in enumerate(lines):
                split_line = line.split()
                if idx == 0:
                    n_atoms = int(line)
                    counter += 1
                elif counter <= n_atoms and counter > 0:
                    elements.append(split_line[0])
                    coordinates.append(
                        [
                            float(split_line[1]),
                            float(split_line[2]),
                            float(split_line[3]),
                        ]
                    )
                    counter += 1
                elif counter == 0:
                    list_of_comments.append(line.replace("\n", ""))
                    counter += 1
                elif counter == n_atoms + 1:
                    counter = -1
                    list_of_coordinates.append(np.array(coordinates))
                    list_of_elements.append(np.array(elements))
                    n_atoms = int(line)
                    coordinates = []
                    elements = []
                    counter = 0
        return list_of_coordinates, list_of_elements, list_of_comments


def main():
    """
    Main function that is called when the script is executed via ensplit."""
    import argparse

    # Step 1: Create the parser
    parser = argparse.ArgumentParser(
        description="Helper script to split an ensemble file into single ensembles with same unique molecular barcode."
    )

    # Step 2: Add arguments
    parser.add_argument("file", type=str, help="Input file")
    parser.add_argument(
        "-T", "--threads", type=int, default=1, help="Number of threads"
    )
    parser.add_argument(
        "-p", "--progress", action="store_true", help="Show progress bar"
    )
    parser.add_argument(
        "-topo", "--topo", action="store_true", help="Only evaluate topology"
    )
    # Step 3: Parse the arguments
    args = parser.parse_args()
    if args.topo:
        topo = True
    else:
        topo = False
    filename = args.file
    parent_dir = os.path.dirname(filename)
    if parent_dir == "":
        parent_dir = "."
    es = EnsembleSplitter(filename=filename)
    es.split_ensemble(threads=args.threads, show_progress_bar=args.progress, topo=topo)
    es.print_ensembles(parent_dir)


if __name__ == "__main__":
    main()
