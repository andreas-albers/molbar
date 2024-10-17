import sys
import os
import json
from molbar.helper.debug import _get_debug_path, _mkdir_debug_path
from molbar.helper.parser import MolParser
from molbar.helper.printer import Printer
from molbar.idealize import idealize_structure_from_file
from molbar.barcode import get_molbar_from_file, get_molbars_from_files


def main() -> None:
    """
    Main function for molbar. Parses arguments and calls the appropriate functions.
    """
    arguments = MolParser().return_arguments()
    if arguments["mode"] == "opt":
        _init_idealization(arguments)
    else:
        _init_barcode_calculation(arguments)


def _init_idealization(arguments: dict) -> None:
    """
    Initializes the idealization of a structure from a file.
    Args:
        arguments (dict): A dictionary containing the parsed arguments from the commandline.
    Raises:
        NotImplementedError: If multiple files are provided.
    """
    if len(arguments["files"]) > 1:
        raise NotImplementedError(
            "Idealization of multiple files is not implemented yet."
        )
    result = idealize_structure_from_file(
        arguments["files"][0],
        return_data=arguments["data"],
        timing=arguments["time"],
        input_constraint=arguments["constraints"][0],
        write_trj=True,
    )
    if arguments["data"]:
        debug_path = _get_debug_path(arguments["files"][0])
        with open(os.path.join(debug_path, "output.json"), "w") as f:
            json.dump(result[4], f, indent=4)
    n_atoms = result[3]
    energy = result[0]
    coordinates = result[1]
    elements = result[2]
    filepath = (
        arguments["files"][0]
        .replace(".xyz", "")
        .replace(".mol", "")
        .replace(".sdf", "")
        .replace(".coord", "")
        + ".opt"
    )
    Printer(
        n_atoms=n_atoms,
        energy=energy,
        coordinates=coordinates,
        elements=elements,
        path=filepath,
    ).print()


def _init_barcode_calculation(arguments: dict) -> None:
    """
    Initializes the calculation of the molecular barcode.
    Args:
        arguments (dict): A dictionary containing the parsed arguments from the commandline.
    """
    if len(arguments["files"]) == 1:
        if arguments["data"]:
            result = get_molbar_from_file(
                arguments["files"][0],
                total_charge=arguments["charge"],
                return_data=arguments["data"],
                timing=arguments["time"],
                input_constraint=arguments["constraints"][0],
                mode=arguments["mode"],
                write_trj=True,
            )
            barcode = result[0]
            debug_path = _get_debug_path(arguments["files"][0])
            if not os.path.isdir(debug_path):
                _mkdir_debug_path(debug_path)
            with open(os.path.join(debug_path, "output.json"), "w") as f:
                json.dump(result[1], f, indent=4)
        else:
            barcode = get_molbar_from_file(
                arguments["files"][0],
                total_charge=arguments["charge"],
                return_data=arguments["data"],
                timing=arguments["time"],
                input_constraint=arguments["constraints"][0],
                mode=arguments["mode"],
            )
        if arguments["save"]:
            file_without_extension = os.path.splitext(arguments["files"][0])[0]
            output = file_without_extension + ".mb"
            with open(output, "w") as output_file:
                print(barcode, file=output_file)
        else:
            print(barcode, file=sys.stdout)

    else:

        if arguments["data"]:
            results = get_molbars_from_files(
                arguments["files"],
                total_charges=arguments["charge"],
                return_data=arguments["data"],
                threads=arguments["threads"],
                timing=arguments["time"],
                input_constraints=arguments["constraints"],
                progress=arguments["progress"],
                mode=arguments["mode"],
                write_trj=True,
            )
            barcodes = [result[0] for result in results]
            for filename, result in zip(arguments["files"], results):
                debug_path = _get_debug_path(filename)
                if not os.path.isdir(debug_path):
                    _mkdir_debug_path(debug_path)
                with open(os.path.join(debug_path, "output.json"), "w") as f:
                    json.dump(result[1], f, indent=4)

        else:
            barcodes = get_molbars_from_files(
                arguments["files"],
                total_charges=arguments["charge"],
                return_data=arguments["data"],
                threads=arguments["threads"],
                timing=arguments["time"],
                input_constraints=arguments["constraints"],
                progress=arguments["progress"],
                mode=arguments["mode"],
            )
        if arguments["save"]:
            for filename, barcode in zip(arguments["files"], barcodes):
                file_without_extension = os.path.splitext(filename)[0]
                output = file_without_extension + ".mb"
                with open(output, "w") as output_file:
                    print(barcode, file=output_file)
        else:
            for barcode in barcodes:
                print(barcode, file=sys.stdout)


if __name__ == "__main__":
    main()
