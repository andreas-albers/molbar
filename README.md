
# Molbar - GrebGroup Version

Customized version of [MolBar](https://git.rwth-aachen.de/bannwarthlab/molbar) v1.1.3:

* The files ideal_geometries.json, soap.json, and element_data.json are loaded only once per MolGraph constructor call and not multiple times during MolBar calculation, reducing overhead.
* Custom topology can be passed to `get_molbar_from_coordinates` using the keyword `edges`. This allows to override MolBar's edge detection, which only provides meaningful results for 2c2e bonds.
* Custom bond orders can be passed to `get_molbar_from_coordinates` using the keyword `bond_orders`. This allows to override MolBar's bond order detection, which only provides meaningful results for 2c2e bonds.
* Formal charges can be passed to `get_molbar_from_coordinates` using the keyword `formal_charges`, enabeling a more precise VSEPR detection by considering lone pairs.
* If the custom VSEPR detection is enabled, the keyword `stereo_settings` allows controlling stereochemistry detection. Default settings are:
  ```python
  {'amine_inversion': True,
  'imine_inversion': True,
  'ring_planarisation': True,
  'ring_hetero_planarisation': True}
  ```
  * `amine_inversion`: planarisation of amines -> amine stereochemistery is ignored
  * `imine_inversion`: linearisation of imines -> imine stereochemistry is ignored
  * `ring_planarisation`: ring structures are planarised -> ring sterochemistry is ignored
  * `ring_hetero_planarisation`: ring atoms with one non-ring atom neighbor are always planarised

As the molecular topology is known in some use cases (_e.g._ comparison of confomer ensembles of a known structure), the above modifications allow reliable MolBar generation of inorganic compounds.

```python
 from molbar.barcode import get_molbar_from_coordinates
```
```python
  def get_molbar_from_coordinates(coordinates: list, elements: list, edges=None, bond_orders=None, formal_charges=None, stereo_settings=None, return_data=False, timing=False, input_constraint=None, mode="mb") -> Union[str, dict]
```
Custom Arguments:
```
  edges (list): Molecular topology provided as a connectivity list, e.g.:
  [[0, 1], [0, 2], [2, 3], [2, 4]]
```
```
  bond_orders (dict): A dict with bond orders, e.g.:
  {'single': [[0, 2], [2, 3]], 'double': [[0, 1], [2, 4]], 'triple': []} 
```
```
  formal_charges (list): A list of formal charges, ordered by atom index, e.g.:
  [-1, 0, 0, 0, 0]
```
```
  stereo_settings (dict): A dict for controlling stereochemistry, e.g.:  {'amine_inversion': True, 'imine_inversion': True, 'ring_planarisation': True, 'ring_hetero_planarisation': True}
```
```
Returns:
  Union[str, dict]: Either MolBar or the MolBar and MolBar data.
```

# MolBar

<div align="center">
<img src="logo.png" alt="logo" width="400" />
</div>

This package offers an implementation of the Molecular Barcode (MolBar), a molecular identifier inspired by quantum chemistry, designed to ensure data uniqueness in chemical structure databases. The identifier is optimized for computational chemistry applications, such as black-box chemical space exploration, by generating the identifier directly from 3D Cartesian coordinates and encoding the 3D shape of the molecule, including information typically not found in a Molfile. This encompasses details about highly distorted structures that deviate significantly from VSEPR geometries. It supports both organic and inorganic molecules and aims to describe relative and absolute configurations, including centric, axial/helical, and planar chirality.


[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![PyPI - Version](https://img.shields.io/pypi/v/molbar.svg)](https://pypi.org/project/molbar)
[![Python Versions](https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)](https://pypi.org/project/molbar)
[![Downloads](https://static.pepy.tech/badge/molbar)](https://pepy.tech/project/molbar)
[![Downloads](https://static.pepy.tech/badge/molbar/month)](https://pepy.tech/project/molbar)
[![ChemRxiv Paper](https://img.shields.io/badge/DOI-10.1038%2Fs41586--020--2649--2-blue)](
https://chemrxiv.org/engage/chemrxiv/article-details/65e3cd80e9ebbb4db9c71da0)

- **ChemRxiv Paper:** https://chemrxiv.org/engage/chemrxiv/article-details/65e3cd80e9ebbb4db9c71da0
- **Documentation:** https://git.rwth-aachen.de/bannwarthlab/molbar/-/blob/main/README.md?ref_type=heads
- **Source code:** https://git.rwth-aachen.de/bannwarthlab/molbar
- **Bug reports:** https://git.rwth-aachen.de/bannwarthlab/molbar/-/issues
- **Email contact:** van.staalduinen@pc.rwth-aachen.de

It does this by fragmentating a molecule into rigid parts which are then idealized with a specialized non-physical force field. The molecule is then described by different matrices encoding topology (connectivity), topography (3D positions of atoms after input unification), and absolute configuration (by calculating a chirality index). The final barcode is the concatenated spectra of these matrices.

## Current Limitations 

**The input file must contain 3D coordinates and explicit hydrogens.**

MolBar should work well for organic and inorganic molecules with typical 2c2e bonding. It can describe molecules based on their relative and absolute configuration, including centric, axial/helical and planar chirality.

Given that 3D Cartesian coordinates are the usual starting point, problems may arise when determining which atoms are bonded, particularly in metal complexes with η-bonds. Additionally, challenges can occur if the geometry around a metal in a complex cannot be classified by one of the standard VSEPR models. Solutions to these issues are being developed and are released with the next versions. If you encounter difficulties, you can use the -d option when using MolBar as a command-line tool or set write_trj=True when using MolBar as a Python module to examine the optimized trajectories of each fragment. If anything is unclear or if you encounter any unusual behavior, please report it by posting issues or via email at van.staalduinen@pc.rwth-aachen.de.

For rigidity analysis, MolBar only considers double/triple bonds and rings to be rigid. For example, an obstacle to rotation due to bulkiness of substituents is not taken into account, but can be added manually from the input file (additional dihedral constraint, but that should be used as an exception and carefully).


## Getting started (tested on Linux and macOS, compiling works for Windows only in WSL)

### For Linux/macOS

Using a virtual environment is highly recommended because it allows you to create isolated environments with their own dependencies, without interfering with other Python projects or the system Python installation. This ensures that your Python environment remains consistent and reproducible across different machines and over time. To create one, type in the following command in your terminal:

```bash
 python3 -m venv your_path
```
To activate the environment, type in:
```bash
 source your_path/bin/activate
```
To install Molbar, enter the following command in your terminal:

```bash
pip install molbar
```

### For Windows

Since compiling in a standard Windows environment does not work yet, it is highly recommended to use the WSL (Windows Subsystem for Linux) extension. Simply follow this installation guide: https://learn.microsoft.com/en-us/windows/wsl/install. Note that a Fortran compiler needs to be installed manually in the WSL environment. Otherwise, the installation of MolBar will result in an error.

For Python usage, it is highly recommended to use Visual Studio Code (VSC) as it provides specific extensions to code directly in WSL. A more detailed guide can be found here: https://code.visualstudio.com/docs/remote/wsl


# MolBar Structure

For l-alanine, the MolBar reads:

```text

MolBar | 1.1.3 | C3NO2H7 | 0 | -339 -140 -110 -32 13 20 20 20 160 237 432 528 850 | -209 -8 130 160 354 633 | -108 -79 -42 -24 11 11 11 47 75 140 293 433 891 | -11 0 0 0 31

```

The different parts of MolBar are defined as follows:

```text
Version: 1.1.3
Molecular Formula: C3H7NO2 
Topology Spectrum: -339 -140 -110 -32 13 20 20 20 160 237 432 528 850 (Encoding atomic connectivity)
Heavy Atom Topology Spectrum: -209 -8 130 160 354 633 (Encoding atomic connectivity without hydrogen atoms. So if for two molecules, the topology spectra are different but the tautomer spectra are the same, both molecules are tautomeric structures)
Topography Spectrum : -108 -79 -42 -24 11 11 11 47 75 140 293 433 891 (3D arrangement of atoms in Cartesian space, also describes diastereomerism)
Chirality: -11 0 0 0 31 (Encoding absolute configuration for each fragment)
```

The chirality barcode can only be compared between two molecules if their topology and topography barcodes are identical. If the chirality barcodes differ only in their sign, this indicates that the two molecules are enantiomers.

# Generating MolBar

MolBar can be generated using either of the following methods:

	1.	Python Interface: Refer to the Python Module Usage for detailed instructions on generating MolBar using Python.

	2.	Command Line Interface with the following command: molbar coord.xyz (see Command Line Usage for more information)


## Python Module Usage

MolBar can be generated by Python function calls:

1. for a single molecule with ```get_molbar_from_coordinates``` by specifying the Cartesian coordinates as a list,
2. for several molecules at once with ```get_molbars_from_coordinates``` by giving a list of lists with Cartesian coordinates,
3. for a single molecule with ```get_molbar_from_file``` by specifying a file path,
4. for several molecules at once with ```get_molbars_from_files``` by specifying a list of file paths.

### 1. get_molbar_from_coordinates
NOTE:
If you need to process multiple molecules at once, it is recommended to use ```get_molbars_from_coordinates``` instead of ```get_molbar_from_coordinates```.

```python
  from molbar.barcode import get_molbar_from_coordinates

  def get_molbar_from_coordinates(coordinates: list, elements: list, return_data=False, timing=False, input_constraint=None, mode="mb") -> Union[str, dict]

```
#### Arguments:

```python
  coordinates (list): Molecular geometry provided by atomic Cartesian coordinates with shape (n_atoms, 3).
```
```python
  elements (list):A list of elements in that molecule. Either the element symbols or atomic numbers can be used.
```
```python
  return_data (bool): Whether to return MolBar data. MolBar can return detailed data, including information about bonding, VSEPR geometries, fragments, and more. This data is useful for understanding what MolBar recognizes and can be leveraged for other projects.
```
```python
  timing (bool): Whether to print the duration of this calculation.
```
```python
  input_constraint (dict, optional): A dict of extra constraints for the calculation. See below for more information. USED ONLY IN EXCEPTIONAL CASES.
```
```python
  mode (str): Whether to calculate the molecular barcode ("mb") or only the topology part of the molecular barcode ("topo").
```

#### Returns:

```python
  Union[str, dict]: Either MolBar or the MolBar and MolBar data.
```

Input constraints should be used only in exceptional cases. However, they may be useful to constrain the molecule with additional dihedrals for rotations that are normally considered around single bonds but whose rotation is hindered (e.g., 90° binol systems with bulky substituents).

  ```python
  {
  'constraints': {
                  'dihedrals': [{'atoms': [1,2,3,4], 'value':90.0},...]} #atoms: list of atoms that define the dihedral, value is the ideal dihedral angle in degrees, atom indexing starts with 1.
  }
  ```


### 2. get_molbars_from_coordinates
NOTE:
If you need to process multiple molecules at once, it is recommended to use this function and specify the number of threads that can be used to process multiple molecules simultaneously.

  ```python
  from molbar.barcode import get_molbars_from_coordinates

  def get_molbars_from_coordinates(list_of_coordinates: list, list_of_elements: list, return_data=False, threads=1, timing=False, input_constraints=None, progress=False,  mode="mb") -> Union[list, Union[str, dict]]:
```

#### Arguments:

```python
  list_of_coordinates (list): A list of molecular geometries provided by atomic Cartesian coordinates with shape (n_molecules, n_atoms, 3).
```
```python
  list_of_elements (list): A list of element lists for each molecule in the list_of_coordinates with shape (n_molecules, n_atoms).  Either the element symbols or atomic numbers can be used.
```
```python
  return_data (bool): Whether to return MolBar data. MolBar can return detailed data, including information about bonding, VSEPR geometries, fragments, and more. This data is useful for understanding what MolBar recognizes and can be leveraged for other projects.
```
```python
  threads (int): Number of threads to use for the calculation. If you need to process multiple molecules at once, it is recommended to use this function and specify the number of threads that can be used to process multiple molecules simultaneously.
```
```python
  timing (bool):  Whether to print the duration of this calculation.
```
```python
  input_constraints (list, optional): A list of constraints for the calculation. Each constraint in that list is a Python dict as shown above for get_molbar_from_coordinates.
```
```python
  progress (bool): Whether to show a progress bar.
```
```python
  mode (str): Whether to calculate the molecular barcode ("mb") or the topology part of the molecular barcode ("topo").
```

#### Returns:

```python
  Union[list, Union[str, dict]]: Either MolBar or the MolBar and MolBar data.
```

### 3. get_molbar_from_file
NOTE:
If you need to process multiple molecules provided by files at once, it is recommended to use ```get_molbars_from_files``` instead of ```get_molbar_from_file```.

  ```python
  from molbar.barcode import get_molbar_from_file

  def get_molbar_from_file(file: str, return_data=False, timing=False, input_constraint=None, mode="mb", write_trj=False) -> Union[str, dict]:
```
#### Arguments:

```python
  file (str): The path to the file containing the molecule information. The file can be in formats such as XYZ, Turbomole coord, SDF/MOL, CIF, PDB, Gaussian (.com, .gjf, .log), and many others. For a complete list, refer to the ASE documentation, ASE I/O Formats: https://wiki.fysik.dtu.dk/ase/ase/io/io.html). The file must contain 3D coordinates and explicit hydrogen atoms.
```
```python
  return_data (bool): Whether to return MolBar data. MolBar can return detailed data, including information about bonding, VSEPR geometries, fragments, and more. This data is useful for understanding what MolBar recognizes and can be leveraged for other projects.
```
```python
  timing (bool): Whether to print the duration of this calculation.
```
```python
  input_constraint (dict, optional): A dict of extra constraints for the calculation. See below for more information. USED ONLY IN EXCEPTIONAL CASES.
```
```python
  mode (str): Whether to calculate the molecular barcode ("mb") or only the topology part of the molecular barcode ("topo").
```
```python
  write_trj (bool, optional): Whether to write a trajectory of the unification process. Defaults to False.
```
    
#### Returns:

```python
  Union[str, dict]: Either MolBar or the MolBar and MolBar data.
  ```

Example for input file in .yml format. Input constraint should be used only in exceptional cases. However, it may be useful to constrain bonds with a additional dihedral for the barcode that are normally considered single bonds but whose rotation is hindered (e.g., 90° binol systems with bulky substituents).
```yml
constraints:
  dihedrals:
    - atoms: [30, 18, 14, 13]  # List of atoms involved in the dihedral
      value:  90.0  # Actual values for the dihedral parameters
```


### 4. get_molbars_from_files

NOTE:
If you need to process multiple molecules at once, it is recommended to use this function and specify the number of threads that can be used to process multiple molecules simultaneously.

  ```python
  from molbar.barcode import get_molbars_from_files

  def get_molbars_from_files(files: list, return_data=False, threads=1, timing=False, input_constraints=None, progress=False, mode="mb", write_trj=False) ->Union[list, Union[str, dict]]:
  ```

#### Arguments:
```python
  files (list): The list of paths to the files containing the molecule information. The file can be in formats such as XYZ, Turbomole coord, SDF/MOL, CIF, PDB, Gaussian (.com, .gjf, .log), and many others. For a complete list, refer to the ASE documentation, ASE I/O Formats: https://wiki.fysik.dtu.dk/ase/ase/io/io.html). The file must contain 3D coordinates and explicit hydrogen atoms.
```
```python
  return_data (bool): Whether to return MolBar data. MolBar can return detailed data, including information about bonding, VSEPR geometries, fragments, and more. This data is useful for understanding what MolBar recognizes and can be leveraged for other projects.
```
```python
  threads (int): Number of threads to use for the calculation. If you need to process multiple molecules at once, it is recommended to use this function and specify the number of threads that can be used to process multiple molecules simultaneously.
```
```python
  timing (bool):  Whether to print the duration of this calculation. 
```
```python
  input_constraints (list, optional): A list of file paths to the input files for the calculation. Each constrained is specified by a file path to a .yml file, as shown above for get_molbar_from_file.
```
```python
  progress (bool): Whether to show a progress bar.
```
```python
  mode (str): Whether to calculate the molecular barcode ("mb") or the topology part of the molecular barcode ("topo").
```
```python
  write_trj (bool, optional): Whether to write a trajectory of the unification process. Defaults to False.
```

#### Returns:
```python
      Union[list, Union[str, dict]]: Either MolBar or the MolBar and MolBar data.

  ```


## Commandline Usage

MolBar can also be used as commandline tool. Just simply type:

```
molbar coord.xyz
```
and the MolBar is printed to the stdout.

NOTE:
If you need to process several molecules at once, it is recommended to pass all molecules to the code at once (e.g. with *.xyz) while specifying the number of threads the code should use:
```bash
molbar *.xyz -T N_threads -s
```
The latter option (-s) is used to store the barcode to .mb files. 


Further, the commandline tool provides several options:

```bash
Basic usage example: 
molbar coord.xyz

Usage with multiple files: 
molbar *.xyz -T 8

Usage for topology barcode only: 
molbar coord.mol -m topo

Usage to return additional data: 
molbar coord.pdb -d

Usage to print timings: 
molbar coord.com -t
```
#### Arguments and Options:

```bash
positional arguments:
  file(s)    

  The file can be in formats such as XYZ, Turbomole coord, SDF/MOL, CIF, PDB, Gaussian (.com, .gjf, .log), and many others. For a complete list, refer to the ASE documentation, ASE I/O Formats: https://wiki.fysik.dtu.dk/ase/ase/io/io.html). The file must contain 3D coordinates and explicit hydrogen atoms.
```
```bash
  -m {mb,topo,opt}, --mode {mb,topo,opt}

  The mode to use for the calculations (either "mb" (default, calculates MolBar), "topo" (topology part only)
  or "opt" (using stand-alone force field idealization, writes ".opt" with final structure))
```
```bash
  -i path/to/input --inp path/to/input

  Path to input file in .yml format to add further constraints. Example input can be found below.
```
```bash
  -d, --data           Whether to print MolBar data. 
                        Writes a "filename/" directory containing a json file with
                        important information that defines MolBar. Writes idealization trajectories of each fragment to same directory.
```
```bash
  -T number_of_threads, --threads number_of_threads
                        The number of threads to use for parallel processing of several files. MolBar generation for a single file is not parallelized. Should be used together with -s/--save (e.g. molbar *.xyz -T 8 -s)
```
```bash
  -s, --save            Whether to save the result to a file of type "filename.mb"
```
```bash
  -t, --time            Print out timings.
```
```bash
  -p, --progress        Use a progress bar when several files are handled.
```

Example for input file constraints in yml format. Input constraint should be used only in exceptional cases. However, it may be useful to constrain bonds with a additional dihedral for the barcode that are normally considered single bonds but whose rotation is hindered (e.g., 90° binol systems with bulky substituents).

```yml
constraints:
  dihedrals:
    - atoms: [30, 18, 14, 13]  # List of atoms involved in the dihedral
      value:  90.0  # Actual values for the dihedral parameters
```



## Using the unification force field for the whole molecule.

The force field can be used to idealize the structure of a whole molecule where the coordinates are either given in Python by a file:

1. as a command line tool with the ```molbar coord.xyz -m opt``` option
2. in Python with ```idealize_structure_from_file``` by providing a file path
3. in Python with ```idealize_structure_from_coordinates``` by providing Cartesian coordinates as a list
 

### Commandline tool
```text
molbar coord.xyz -m opt
```
This writes a coord.opt file that contains the idealized coordinates.

### In Python from a file:
```python
  from molbar.idealize import idealize_structure_from_file

  def idealize_structure_from_file(file: str, return_data=False, timing=False, input_constraint=None,  write_trj=False) -> Union[list, str]
```
#### Arguments:
```python
  file (str): The path to the input file to be processed. 
  
  The file can be in formats such as XYZ, Turbomole coord, SDF/MOL, CIF, PDB, Gaussian (.com, .gjf, .log), and many others. For a complete list, refer to the ASE documentation, ASE I/O Formats: https://wiki.fysik.dtu.dk/ase/ase/io/io.html). The file must contain 3D coordinates and explicit hydrogen atoms.
```
```python
  return_data (bool): Whether to return MolBar data. MolBar can return detailed data, including information about bonding, VSEPR geometries, fragments, and more. This data is useful for understanding what MolBar recognizes and can be leveraged for other projects.
```
```python
  timing (bool): Whether to print the duration of this calculation.
```
```python
  input_constraint (str): The path to the input file containing the constraint for the calculation. See down below for more information.
```
```python
  write_trj (bool, optional): Whether to write a trajectory of the unification process. Defaults to False.
```
#### Returns:
```python
  n_atoms (int): Number of atoms in the molecule.
```
```python
  energy (float): Final energy of the molecule after idealization.
```
```python
  coordinates (list): Final coordinates of the molecule after idealization.
```
```python
  elements (list): Elements of the molecule.
```
```python
    data (dict): Molbar data.
```

This is an example input as a yml file:
```yml
bond_order_assignment: False  # False if bond order assignment should be skipped, only reasonable opt mode (standalone force-field optimization)
cycle_detection: True # False if cycle detection should be skipped, only reasonable opt mode (standalone force-field optimization).
repulsion_charge: 100.0 # Charge used for the Coulomb term in the Force field, every atom-atom pair uses the same charge, only reasonable opt mode (standalone force-field optimization). Defaults to 100.0
set_edges: True #False if no bonds should be constrained automatically.
set_angles: True #False if no angles should be constrained automatically.
set_dihedrals: True # False if no dihedrals should be constrained automatically.
set_repulsion: True #False if no coulomb term should be used automatically.

constraints:
  bonds:
    - atoms: [19, 23]  # List of atoms involved in the bond
      value: 1.5  # Ideal bond length. 
  angles:
    - atoms: [19, 23, 35]  # List of atoms involved in the angle
      value: 45.0  # Angle to which the angle between the three atoms is to be constrained
    - atoms: [35, 23, 19]  # List of atoms involved in the angle
      value: 45.0  # Angle to which the angle between the three atoms is to be constrained

  dihedrals:
    - atoms: [30, 18, 14, 13]  # List of atoms involved in the dihedral
      value:  90.0  # Actual values for the dihedral parameters
```

### In Python from a list of Cartesian coordinates:
```python
from molbar.idealize import idealize_structure_from_coordinates

def idealize_structure_from_coordinates(coordinates: list, elements: list, return_data=False, timing=False, input_constraint=None) -> Union[list, str]:
```
#### Arguments:

```python
  coordinates (list): Cartesian coordinates of the molecule.
```
```python
  elements (list): Elements of the molecule. Either the element symbols or atomic numbers can be used.
```
```python
  return_data (bool, optional): Whether to return MolBar data. MolBar can return detailed data, including information about bonding, VSEPR geometries, fragments, and more. This data is useful for understanding what MolBar recognizes and can be leveraged for other projects.
```
```python
  timing (bool, optional): Whether to print the duration of this calculation.
```
```python
  input_constraint (dict, optional): The constraint for the calculation. See documentation for more information.
```
  
#### Returns:
```python
  n_atoms (int): Number of atoms in the molecule.
```
```python
  energy (float): Final energy of the molecule after idealization.
```
```python
  coordinates (list): Final coordinates of the molecule after idealization.
```
```python
  elements (list): Elements of the molecule.
```
```python
    data (dict): Molbar data.
```

This is an example input as a Python dict:

```python
  {'bond_order_assignment': True, #False if bond order assignment should be skipped, only reasonable opt mode (standalone force-field optimization)
  'cycle_detection': True, #False if cycle detection should be skipped, only reasonable opt mode (standalone force-field optimization).
  'set_edges': True #False if no bonds should be constrained automatically.
  'set_angles': True #False if no angles should be constrained automatically.
  'set_dihedrals': True #False if no dihedrals should be constrained automatically.
  'set_repulsion': True #False if no coulomb term should be used automatically.
  'repulsion_charge': 100.0, # Charged used for the Coulomb term in the Force field, every atom-atom pair uses the same charge, only reasonable opt mode (standalone force-field optimization). Defaults to 100.0
  'constraints': {'bonds': [{'atoms': [1,2], 'value':1.5},...], #atoms: list of atoms that define the bond, value is the ideal bond length in angstrom, atom indexing starts with 1.
                  'angles': [{'atoms': [1,2,3], 'value':90.0},...], #atoms: list of atoms that define the angle, value is the ideal angle in degrees, atom indexing starts with 1.
                  'dihedrals': [{'atoms': [1,2,3,4], 'value':180.0},...]}  #atoms: list of atoms that define the dihedral, value is the ideal dihedral angle in degrees, atom indexing starts with 1.
  }
```

# Additional Command line Scripts provided by the MolBar package.

## Ensplit
Helper script to split an ensemble file into individual ensembles, each with a unique molecular barcode:

```bash
ensplit ensemble.trj -T 8 -topo
```

```bash
file                  Input file 

The file can be in formats such as XYZ, Turbomole coord, SDF/MOL, CIF, PDB, Gaussian (.com, .gjf, .log), and many others. For a complete list, refer to the ASE documentation, ASE I/O Formats: https://wiki.fysik.dtu.dk/ase/ase/io/io.html). The file must contain 3D coordinates and explicit hydrogen atoms.

```
```bash
  -T THREADS, --threads number_of_threads. Number of threads for processing multiple files.
```
```bash
  -p, --progress        Show progress bar
```
```bash
  -topo, --topo         Only evaluate topology
```

## Princax

Helper script to align the molecule to the principal axes.

```bash
princax coord.xyz -r
```

```bash
file                  Input file 

The file can be in formats such as XYZ, Turbomole coord, SDF/MOL, CIF, PDB, Gaussian (.com, .gjf, .log), and many others. For a complete list, refer to the ASE documentation, ASE I/O Formats: https://wiki.fysik.dtu.dk/ase/ase/io/io.html). The file must contain 3D coordinates and explicit hydrogen atoms.

```
```bash
-r, --replace  Overwrite the input file with the aligned coordinates, or print to a new file if not specified.
```

## Invstruc

Helper script to invert structures to yield enantiomers.

```bash
invstruc coord.xyz
```

```bash
file                  Input file 

The file can be in formats such as XYZ, Turbomole coord, SDF/MOL, CIF, PDB, Gaussian (.com, .gjf, .log), and many others. For a complete list, refer to the ASE documentation, ASE I/O Formats: https://wiki.fysik.dtu.dk/ase/ase/io/io.html). The file must contain 3D coordinates and explicit hydrogen atoms.

```

## Acknowledgments

MolBar relies on the following libraries
and packages:

*   [ase](https://wiki.fysik.dtu.dk/ase/)
*   [dscribe](https://singroup.github.io/dscribe/latest/)
*   [networkx](https://networkx.org/)
*   [NumPy](https://numpy.org)
*   [Numba](https://numba.pydata.org)
*   [SciPy](https://scipy.org)
*   [tqdm](https://github.com/tqdm/tqdm)
*   [joblib](https://joblib.readthedocs.io/en/latest/)

Thank you!


## License and Disclaimer

MIT License

Copyright (c) 2022 Nils van Staalduinen, Christoph Bannwarth

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
