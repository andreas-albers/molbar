# Changelog

## [1.1.2] - 2024-07-19

MolBar 1.1.2 is a maintenance release that fixes bugs discovered after the 1.1.1 release and focuses mainly on robustness and stability.

### API Change:
- **Breaking Change**: Moved `idealize_structure_from_file` and `idealize_structure_from_coordinates` from `molbar.barcode` to `molbar.idealize`.
- **Additive Change**: Instead of element symbols, atomic numbers can now be provided to the elements parameter in `get_molbar_from_coordinates` and `get_molbars_from_coordinates`. So ```z = [6, 1, 1, 1, 1] ``` instead of ```["C", "H", "H", "H", "H"]``` is also valid. See documentation for more details.
- **Additive Change**: The `molbar` package now includes three scripts that can be used as command-line tools. The latter two scripts are commonly used during the development of MolBar, which is why they have been included:
  1. `ensplit ensemble.xyz` - Splits an ensemble file into multiple files, each containing structures with the same unique MolBar.
  2. `princax coord.xyz` - Aligns the structure to its principal axes.
  3. `invstruc coord.xyz` - Mirrors the structure to yield the enantiomer.

### Main Changes:
- **Added Total Charge**: Introduced an additional charge layer to MolBar, which represents the total charge of the molecule and is located after the molecular formula. For example, methane cation is represented as `MolBar | 1.1.2 | CH3 | 1 | -121 20 20 381 | 60 | -55 12 12 332 | 0`, while methane anion is represented as `MolBar | 1.1.2 | CH3 | -1 | -121 20 20 381 | 60 | -55 12 12 332 | 0`. See documentation for more details.
- **Resolved Issue #152**: Enhanced error management and type annotations. All inputs are now validated for correctness.
- Ensure wheel compatibility across different platforms using `auditwheel` for Linux and `delocate` for macOS. This includes bundling the necessary Fortran libraries to ensure that all dependencies are included in the wheel. `auditwheel` repairs the wheels by including external shared libraries, copying these libraries into the wheel itself, and automatically modifying the appropriate RPATH entries so that these libraries will be picked up at runtime. Similarly, `delocate` performs these tasks for macOS, ensuring the wheel is self-contained and all dynamic libraries are correctly referenced.
- Added wheels for macOS 11.0+ for x86_64 architecture.
- Enhanced numerical stability for the chirality index by improving the determination of asymmetry in molecules through symmetry operations recognition.
- Improved code quality by performing a comprehensive code cleanup and starting to increase the test suite coverage.
- Adopted the use of the Black linter for consistent code formatting.


## [1.1.1] - 2024-03-14

MolBar 1.1.1 is a maintenance release that fixes bugs discovered after the 1.1.0 release.

Main bug fixes:

- Includes a richer setup.py, since cp3.8 relies on setup.py for package information, as opposed to cp3.9 and newer, which prioritizes information from pyproject.toml.

## [1.1.0] - 2024-03-10
The MolBar v1.1.0 release adds Python 3.12.0 support, continuing from the v1.0.0 release. It's now easier to install with the standard Python way via PyPi with just "pip install molbar" due to its availability on PyPi.org (https://pypi.org/project/molbar/), i.e. without gitlab in between.
Here are the key updates:

### PyPi Upload and Build Wheels:

MolBar is on PyPi (pip install molbar), using build wheels for quicker installations by avoiding compilation of C or Fortran extensions when possible.

Wheels are a package format for Python: Designed to simplify the distribution and installation of Python software, they replace the traditional egg format.
Faster installation: Wheels allow for precompilation of code, eliminating the need for users to compile the software on their own machines. This results in significantly faster installation times compared to installing from source.
Consistency and compatibility: By including compiled binaries, wheels ensure that software runs consistently across environments and reduce the potential for errors associated with compilation on different system.

### Switch to Scikit-build:

With Python 3.12 dropping numpy.distutils, MolBar moved to the scikit-build/CMake to compile Fortran extensions, which supports more compilers and provides better cross-compilation and dependency management.

### More Input File Types:

MolBar now supports additional file types, including Gaussian (.gjf/.com), Turbomole, CIF, MOL, SDF, and XYZ.
The Python versions supported by this release are 3.8-3.12.