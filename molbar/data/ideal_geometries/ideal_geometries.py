from molbar.io.filereader import FileReader
from molbar.helper.vector import Vector
from molbar.helper.printer import Printer
from os import listdir
import json
import numpy as np

if __name__ == "__main__":

    files = [file for file in listdir(".") if file.endswith(".xyz")]

    ideal_geometries = {}

    for file in files:

        n_atoms, geometry, atoms = FileReader(file).read_file()

        geometry = np.array([x - geometry[0] for x in geometry])

        geometry = np.array(
            [
                x / np.linalg.norm(x) if i != 0 else [0.0, 0.0, 0.0]
                for i, x in enumerate(geometry)
            ]
        )

        print(file)

        print(geometry)

        Printer(n_atoms, 0.0, geometry, atoms, file).print()

        angles = []

        for i in range(1, n_atoms):

            for j in range(1, n_atoms):

                if j > i:

                    nodes = [i, 0, j]

                    alpha = float(round(Vector(geometry[nodes]).calculate_angle()))

                    angles.append({"nodes": nodes, "angle": alpha})

        dihedrals = []

        for i in range(1, n_atoms):

            for j in range(1, n_atoms):

                for k in range(1, n_atoms):

                    if (i > k) and (j != i) and (j != k):

                        nodes = [i, 0, j, k]

                        try:

                            theta = float(
                                round(Vector(geometry[nodes]).calculate_dihedral())
                            )

                        except ValueError:

                            continue

                        if theta == -0.0:

                            dihedral = 0.0

                        dihedrals.append({"nodes": nodes, "dihedral": theta})

        ideal_geometries[file.replace(".xyz", "")] = {
            "angles": angles,
            "dihedrals": dihedrals,
        }

    with open("../ideal_geometries.json", "w") as f:

        json.dump(ideal_geometries, f, indent=4)
