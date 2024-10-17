from ase import Atoms
from molbar.io.filereader import FileReader
from dscribe.descriptors import SOAP
import os
import json


def main():
    cur_dir = "./ideal_geometries/"
    onlyfiles = [f for f in os.listdir(cur_dir) if f.endswith(".xyz")]
    data = {}

    for file in onlyfiles:

        print(file)

        n_atoms, geometry, atoms = FileReader(os.path.join(cur_dir, file)).read_file()

        soap = SOAP(
            species=["C", "H"],
            periodic=False,
            r_cut=3.0,
            n_max=2,
            l_max=1,
        )

        atoms = Atoms("".join(atoms), positions=geometry)

        data[file.replace(".xyz", "")] = soap.create(atoms, centers=[0]).tolist()[0]

    print(data)
    with open("soap.json", "w") as output:

        json.dump(data, output, indent=4)


if __name__ == "__main__":

    main()
