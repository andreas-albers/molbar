from argparse import ArgumentParser


class MolParser(ArgumentParser):
    """
    ArgumentParser-derived objects to collect the arguments given by the commandline.

    Methods:
        return_arguments -> dict: Returns a dict of all given arguments, containing following keys:

                files: list, paths to .xyz-files.
                debug: bool, if debug file should be generated.
                threads: int, number of threads (if several .xyz should be handled in parallel
                omb: bool, if MolBar should be stored in .mb-file
                time: bool, if timings should be printed out.
                progress: bool, use a progress bar to illustrate the progress.
    """

    def __init__(self) -> None:
        """
        Constructs MolParser object.

        Args:
            None
        """

        super().__init__()

        self.add_argument("files", type=str, nargs="+", help="file(s)")

        self.add_argument(
            "-i",
            "--inp",
            type=str,
            nargs="+",
            help="Path to input file to add further constraints. Example input can be found in the documentation.",
        )

        self.add_argument(
            "-d",
            "--data",
            action="store_true",
            default=False,
            help='Whether to write MolBar data. Writes a "filename/" directory containing a json file with important information that defines MolBar. Writes idealization trajectories of each fragment to same directory. ',
        )

        self.add_argument(
            "-T",
            "--threads",
            type=int,
            default=1,
            help="The number of threads to use for parallel processing of several files. MolBar generation for a single file is not parallelized. Usage as molbar *.xyz -T 8 -s",
        )

        self.add_argument(
            "-s",
            "--save",
            action="store_true",
            default=False,
            help='Whether to save the result to a file of type "filename.mb"',
        )

        self.add_argument(
            "-t",
            "--time",
            action="store_true",
            default=False,
            help="Print out timings.",
        )

        self.add_argument(
            "-p",
            "--progress",
            action="store_true",
            default=False,
            help="Use a progress bar to illustrate the progress.",
        )

        self.add_argument(
            "-c",
            "--charge",
            type=int,
            default=0,
            help="The total charge of the molecule(s).",
        )

        self.add_argument(
            "-m",
            "--mode",
            type=str,
            choices=["mb", "topo", "opt"],
            default="mb",
            help='The mode to use for the calculations (either "mb" (default, calculates MolBar), "topo" (topology part only) or "opt" (using stand-aloneforce field idealization, writes ".log" with trajectory and ".opt" with final structure)',
        )

        self.arguments = self.parse_args()

    def return_arguments(self) -> dict:
        """
        return_arguments -> dict: Returns a dict of all given arguments, containing following keys:

                files: list, paths to .xyz-files.
                debug: bool, if debug file should be generated.
                threads: int, number of threads (if several .xyz should be handled in parallel)
                omb: bool, if MolBar should be stored in .mb-file.
                time: bool, if timings should be printed out.
        """

        if self.arguments.inp:

            if len(self.arguments.files) != len(self.arguments.inp):

                raise ValueError("Number of files and control files do not match.")

            constraints = self.arguments.inp

        else:

            empty_constraints = None

            constraints = [empty_constraints] * len(self.arguments.files)

        return {
            "files": self.arguments.files,
            "charge": self.arguments.charge,
            "data": self.arguments.data,
            "threads": self.arguments.threads,
            "save": self.arguments.save,
            "time": self.arguments.time,
            "constraints": constraints,
            "progress": self.arguments.progress,
            "mode": self.arguments.mode,
        }
