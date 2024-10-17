import os
from skbuild import setup

# I want my readme to be part of the setup, so let's read it.

lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = lib_folder + '/requirements.txt'

if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requirements = list(f.read().splitlines())

with open(os.path.join(lib_folder, 'README.md'), 'r') as fh:
    long_description = fh.read()


setup(
    name='molbar',
    author='Nils van Staalduinen',
    author_email='van.staalduinen@pc.rwth-aachen.de',
    version='1.1.3',
    description='Molecular Barcode (MolBar): Molecular Identifier for Organic and Inorganic Molecules',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    url='https://git.rwth-aachen.de/bannwarthlab/MolBar',
    install_requires=install_requirements,
    keywords='graph representation, permuation-invariant, molecular barcode',
    packages=['molbar'],
    cmake_args=['-DSKBUILD=ON'],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License"
    ],
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'molbar = molbar.main:main',
            'ensplit = molbar.helper.ensemble_splitter:main',
            'princax = molbar.helper.symcheck:main',
            "invstruc = molbar.helper.invstruc:main",
            "fragcount = molbar.helper.fragcount:main"]
    }
)