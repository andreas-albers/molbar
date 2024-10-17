import sys


class MolBarGenerationError(Exception):
    """
    MolBarGenerationError, Exception, when Molbar could not be generated for specific file.

    Args:
        file(str): path to file
    """

    def __init__(self, file, *args):

        self.file = file
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "MolBarGenerationError, {0} ".format(self.message)
        else:
            return "MolBarGenerationError, MolBar could not be generated for molecule {0}.".format(
                self.file
            )


class FileNotFound(Exception):
    """
    FileNotFound, Exception, input file was not found.

    Args:
        file(str): path to file
    """

    def __init__(self, file, *args):

        self.file = file
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "FileNotFound, {0} ".format(self.message)
        else:
            return "FileNotFound, File {0} does not exist.".format(self.file)


class FileFormatNotSupported(Exception):
    """
    FileFormatNotSupported, Exception, when Molbar could not be generated for specific file type. Must be either .xyz, .coord/.tmol or .sdf/.mol (for the latter two: it must be 3D coordinates).
    Args:
        file(str): path to file
    """

    def __init__(self, file, *args):

        self.file = file
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "FileFormatNotSupported, {0} ".format(self.message)
        else:
            return "FileFormatNotSupported, MolBar could not be generated for file {0} as it is not a supported file format! The input file is not supported. Please specify the file type with the extensions .xyz, .coord/.tmol and .sdf/.mol (for the latter two: it must be 3D coordinates).".format(
                self.file
            )


class NotXYZFormat(Exception):
    """
    NotXYZFormat, Exception, when Molbar could not be generated for specific file due to input file format not as .xyz.

    Args:
        file(str): path to file
    """

    def __init__(self, file, *args):

        self.file = file
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "NotXYZFormat, {0} ".format(self.message)
        else:
            return "NotXYZFormat, MolBar could not be generated for file {0} as it is not a .xyz!".format(
                self.file
            )


class NotMolFormat(Exception):
    """
    NotMolFormat, Exception, when Molbar could not be generated for specific file due to input file format not as .mol.

    Args:
        file(str): path to file
    """

    def __init__(self, file, *args):

        self.file = file
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "NotMolFormat, {0}".format(self.message)
        else:
            return "NotMolFormat, Molbar could not be generated for file {0} as it is not a .mol!".format(
                self.file
            )


class NotV2000orV3000Format(Exception):
    """
    NotV2000orV3000Format, Exception, when MolBar could not be genrated for specific file due to .mol input file format not as V2000 or as V3000.

    Args:
        file(str): path to file

    """

    def __init__(self, file, *args):

        self.file = file
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "NotV2000Format, {0}".format(self.message)
        else:
            return "NotV2000Format, MolBar could not be generated for file {0} as it is not a .mol file in the V2000 format!".format(
                self.file
            )


class NotCoordFormat(Exception):
    """
    NotCoordFormat, Exception, when Molbar could not be generated for specific file due to input file format not as .coord.

    Args:
        file(str): path to file
    """

    def __init__(self, file, *args):

        self.file = file
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "NotCoordFormat, {0} ".format(self.message)
        else:
            return "NotCoordFormat, MolBar could not be generated for file {0} as it is not a .coord!".format(
                self.file
            )


class NotMBFormat(Exception):
    """
    NotMBFormat, Exception, when Molbar could not be read from given sfile.

    Args:
        file(str): path to file
    """

    def __init__(self, file, *args):

        self.file = file
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "NotMBFormat, {0} ".format(self.message)
        else:
            return "NotMBFormat, MolBar could not be read from file {0} as it is not an .mb-file!".format(
                self.file
            )


class NotInputFormat(Exception):
    """
    NotYMLFormat, Exception, when input file does not contain the right keys. Only "bond_order_assignment", "cycle_detection", "constraints" are allowed.

    Args:
        file(str): path to file
    """

    def __init__(self, file, *args):

        self.file = file
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "NotInputFormat, {0} ".format(self.message)
        else:
            return 'NotInputFormat, input file {0} is not in right format. Input file might not contain the right keys. Only "bond_order_assignment", "cycle_detection", "constraints" are allowed.'.format(
                self.file
            )


class UnusualValence(Exception):
    """
    UnusualValence, Exception, when unsual valence is used for certain atom.

    Args:
        i(int): atomic index
        atom(str): element
        valence(int): valence of ith atom
    """

    def __init__(self, path, i: int, atom: str, valence: int, *args):

        self.path = path
        self._i = i + 1
        self._atom = atom
        self._valence = valence

    def __str__(self):

        return f"Warning for molecule {self.path}: Unsual valence excepted for atom {self._i} ({self._atom}): valence of {self._valence} "


class AtomNotAdjacent(Exception):
    """
    AtomNotAdjacent, Exception, when jth atom determining direction of path in ring is not adjacent to ith atom in cis/trans determination.

    Args:
        i(int): atomic index
        atom(str): element
        valence(int): valence of ith atom
    """

    def __init__(self, i: int, j: int, ring: list, *args):

        self._i = i + 1
        self._j = j + 1
        self._ring = [str(x + 1) for x in ring]

    def __str__(self):

        return f'Warning for molecule {sys.argv[1]}: Atom {self._j} is not adjacent to {self._i} in ring ({" ".join(self._ring)}): Might cause a false cis/trans determination for atom {self._i}.'


class InitialMinimum(Exception):
    """
    NotMBFormat, Exception, when initial structure is minimum.

    Args:
        file(str): path to file
    """

    def __init__(self):

        self.type = "initial"

    def __str__(self):

        return f"InitialMinimum, initial structure is minimum."


class OptConvergence(Exception):
    """
    OptConvergence, Exception, when optimization has converged.

    """

    def __init__(self, type):

        self.type = type

    def __str__(self):

        return f"OptConvergence of type {self.type}"


class NextStep(Exception):
    """
    NextStep, Exception, when next opt step should be carried out.
    """

    def __init__(self, *args):

        if args:

            self.type = args[0]
        else:

            self.type = ""

    def __str__(self):

        return f"NextStep in optimization should follow."


class NoFragments(Exception):
    """
    NoFragments, Exception, when no fragments were detected..
    """

    def __init__(self, *args):

        if args:

            self.type = args[0]
        else:

            self.type = ""

    def __str__(self):

        return f"No fragments were detected."
