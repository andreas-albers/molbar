#!/bin/bash

# Print Python version
echo "Using Python version: $(python --version)"

# Define paths and version variables
#molbar_path=~/molbar
molbar_path=~/molbar_project/source
molbar_version=1.1.2
dist_path=$molbar_path/dist
base_venv_name="venv_3"
path_to_venv=$molbar_path/build_wheels
#homebrew_path=/usr/local/bin
homebrew_path=/opt/homebrew/bin
build_macos_target=14.0
final_macos_target=11.0
#architecture="x86_64"
architecture="arm64"
# Create distribution directory if it doesn't exist
mkdir -p $dist_path

# Deactivate any active virtual environment
if [ -n "$VIRTUAL_ENV" ]; then
    deactivate
fi

# Loop through the Python versions 8 through 12
for version in {8..12}
do
    venv_name="${path_to_venv}/${base_venv_name}${version}"
    
    echo "Creating virtual environment: ${venv_name}"
    
    # Remove existing virtual environment if it exists
    if [ -d "${venv_name}" ]; then
        rm -r "${venv_name}"
    else
        echo "Directory ${venv_name} does not exist."
    fi
    
    # Create a new virtual environment
    ${homebrew_path}/python3.${version} -m venv "${venv_name}"
    
    echo "Activating virtual environment: ${venv_name}"
    source "${venv_name}/bin/activate"
    
    echo "Python version in virtual environment: $(python --version)"
    
    echo "Installing necessary packages"
    pip install --upgrade pip
    pip install numpy scikit-build delocate
    
    # Remove previous build artifacts if they exist
    [ -d $molbar_path/_skbuild ] && rm -r $molbar_path/_skbuild || echo "_skbuild directory does not exist."
    [ -d $molbar_path/molbar.egg-info ] && rm -r $molbar_path/molbar.egg-info || echo "molbar.egg-info directory does not exist."
    
    pip uninstall molbar --yes
    
    # Change to the molbar project directory
    cd $molbar_path
    
    # Set environment variables for building with the macOS 11.0 SDK
    export MACOSX_DEPLOYMENT_TARGET=$build_macos_target

    echo "Building wheel"
    python3 setup.py bdist_wheel

    # Construct wheel file name
    wheelname_old="molbar-${molbar_version}-cp3${version}-cp3${version}-macosx_${build_macos_target%.*}_0_${architecture}.whl"
    wheelname="molbar-${molbar_version}-cp3${version}-cp3${version}-macosx_${final_macos_target%.*}_0_${architecture}.whl"

    # Repair wheel with delocate
    echo "Repairing wheel: ${dist_path}/${wheelname}"
    
    delocate-wheel -v $dist_path/$wheelname_old

    mv $dist_path/$wheelname_old $dist_path/$wheelname

    echo "Testing wheel: ${wheelname}"

    export MACOSX_DEPLOYMENT_TARGET=$final_macos_target
    
    # Install and test the wheel
    pip install $dist_path/$wheelname
    pip install pytest
    pytest .
    
    # Deactivate the virtual environment
    deactivate
    
done

echo "Wheel building and testing completed."