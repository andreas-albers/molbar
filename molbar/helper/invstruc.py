import os
from molbar.molgraph.molgraph import MolGraph
from molbar.helper.printer import Printer


def main():
    """
    Main function that is called when the script is executed via invstruc."""
    import argparse

    # Step 1: Create the parser
    parser = argparse.ArgumentParser(
        description="Helper script to invert structures to yield enantiomers."
    )

    # Step 2: Add arguments
    parser.add_argument("file", type=str, help="Input file")
    # Step 3: Parse the arguments
    args = parser.parse_args()
    filename = args.file
    dirname = os.path.dirname(filename)
    basename = os.path.basename(filename)
    basename_without_extension = os.path.splitext(basename)[0]
    new_filename = os.path.join(dirname, basename_without_extension + "_inverted.xyz")
    molGraph = MolGraph()
    molGraph.from_file(filename)
    coordinates = molGraph.return_node_data(attribute="coordinates")
    elements = molGraph.return_node_data(attribute="elements")
    coordinates[:, 0] *= -1
    Printer(
        n_atoms=len(elements),
        energy=0.0,
        coordinates=coordinates,
        elements=elements,
        path=new_filename,
    ).print()


if __name__ == "__main__":
    main()
