#!/bin/bash

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

echo "Checking and installing dependencies..."

# Update package list
sudo apt update

# Check and install make
if ! command_exists make; then
    echo "Installing make..."
    sudo apt install -y make
else
    echo "make is already installed."
fi

# Check and install GCC
if ! command_exists gcc; then
    echo "Installing GCC..."
    sudo apt install -y gcc
else
    echo "GCC is already installed."
fi

# Check and install OpenMP support
if ! command_exists libgomp1; then
    echo "Installing OpenMP support..."
    sudo apt install -y libgomp1
else
    echo "OpenMP support is already installed."
fi

# Check and install LAPACK and BLAS
if ! ldconfig -p | grep -q liblapack; then
    echo "Installing LAPACK and BLAS..."
    sudo apt install -y liblapack-dev libblas-dev
else
    echo "LAPACK and BLAS are already installed."
fi

# Check and install ARPACK
if ! ldconfig -p | grep -q libarpack; then
    echo "Installing ARPACK..."
    sudo apt install -y libarpack2-dev
else
    echo "ARPACK is already installed."
fi

# Check and install gfortran
if ! command_exists gfortran; then
    echo "Installing gfortran..."
    sudo apt install -y gfortran
else
    echo "gfortran is already installed."
fi

echo "All dependencies are installed and up-to-date."
