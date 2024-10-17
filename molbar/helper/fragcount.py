import networkx as nx
import numpy as np
from molbar.barcode import get_molbar_from_file

def main():
    """
    Main function that is called when the script is executed via fragcount."""
    import argparse

    # Step 1: Create the parser
    parser = argparse.ArgumentParser(
        description="Helper script to count number of molecular fragments."
    )
    # Step 2: Add arguments
    parser.add_argument("file", type=str, help="Input file")
    # Step 3: Parse the arguments
    args = parser.parse_args()
    molbar, data = get_molbar_from_file(args.file, return_data=True,mode="topo")
    n_atoms = len(data["elements"])
    if n_atoms == 0:
        print("There are not atoms in the file.")
        return 
    cn_matrix = np.zeros((n_atoms,n_atoms))
    for bond in data["single_bonds"]: 
        cn_matrix[bond[0]-1, bond[1]-1] = 1.
        cn_matrix[bond[1]-1, bond[0]-1] = 1.
    for bond in data["double_bonds"]: 
        cn_matrix[bond[0]-1, bond[1]-1] = 1.
        cn_matrix[bond[1]-1, bond[0]-1] = 1.
    for bond in data["triple_bonds"]: 
        cn_matrix[bond[0]-1, bond[1]-1] = 1.
        cn_matrix[bond[1]-1, bond[0]-1] = 1.
    G = nx.Graph(cn_matrix)
    n_molecules = len(list(nx.connected_components(G)))
    print(n_molecules)

if __name__ == "__main__":
    main()

