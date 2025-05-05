#!/bin/bash

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to detect package manager
detect_package_manager() {
    if command_exists apt; then
        echo "apt"
    elif command_exists yum; then
        echo "yum"
    elif command_exists pacman; then
        echo "pacman"
    else
        echo "unsupported"
    fi
}

PACKAGE_MANAGER=$(detect_package_manager)

if [ "$PACKAGE_MANAGER" = "unsupported" ]; then
    echo "Unsupported Linux distribution. Exiting."
    exit 1
fi

echo "Detected package manager: $PACKAGE_MANAGER"
echo "Checking and installing dependencies..."

# Update package list
if [ "$PACKAGE_MANAGER" = "apt" ]; then
    sudo apt update
elif [ "$PACKAGE_MANAGER" = "yum" ]; then
    sudo yum makecache
elif [ "$PACKAGE_MANAGER" = "pacman" ]; then
    echo "Refreshing and ranking mirrors using reflector..."
    if ! command_exists reflector; then
        sudo pacman -S --noconfirm reflector
    fi
    sudo reflector --country Poland,Sweden --latest 10 --sort rate --save /etc/pacman.d/mirrorlist
    sudo pacman -Sy
fi

# Function to install a package
install_package() {
    PACKAGE=$1
    if [ "$PACKAGE_MANAGER" = "apt" ]; then
        sudo apt install -y "$PACKAGE"
    elif [ "$PACKAGE_MANAGER" = "yum" ]; then
        sudo yum install -y "$PACKAGE"
    elif [ "$PACKAGE_MANAGER" = "pacman" ]; then
        sudo pacman -S --noconfirm "$PACKAGE"
    fi
}

# Check and install make
if ! command_exists make; then
    echo "Installing make..."
    install_package make
else
    echo "make is already installed."
fi

# Check and install GCC
if ! command_exists gcc; then
    echo "Installing GCC..."
    install_package gcc
else
    echo "GCC is already installed."
fi

# Check and install LAPACK and BLAS
if ! ldconfig -p | grep -q liblapack; then
    echo "Installing LAPACK and BLAS..."
    install_package liblapack-dev
    install_package libblas-dev
else
    echo "LAPACK and BLAS are already installed."
fi

# Check and install ARPACK
if ! ldconfig -p | grep -q libarpack; then
    echo "Installing ARPACK..."
    if [ "$PACKAGE_MANAGER" = "pacman" ]; then
        install_package arpack
    else
        install_package libarpack2-dev
    fi
else
    echo "ARPACK is already installed."
fi

# Check and install gfortran
if ! command_exists gfortran; then
    echo "Installing gfortran..."
    if [ "$PACKAGE_MANAGER" = "pacman" ]; then
        install_package gcc-fortran
    else
        install_package gfortran
    fi
else
    echo "gfortran is already installed."
fi

# Check and install OpenBLAS
if ! ldconfig -p | grep -q libopenblas; then
    echo "Installing OpenBLAS..."
    if [ "$PACKAGE_MANAGER" = "pacman" ]; then
        install_package openblas
    else
        install_package libopenblas-dev
    fi
else
    echo "OpenBLAS is already installed."
fi

# Check and install time
if ! command_exists time; then
    echo "Installing time..."
    if [ "$PACKAGE_MANAGER" = "pacman" ]; then
        install_package time
    elif [ "$PACKAGE_MANAGER" = "apt" ]; then
        install_package time
    elif [ "$PACKAGE_MANAGER" = "yum" ]; then
        install_package time
    fi
else
    echo "time is already installed."
fi

# Check and install valgrind
if ! command_exists valgrind; then
    echo "Installing valgrind..."
    install_package valgrind
else
    echo "valgrind is already installed."
fi

echo "All dependencies are installed and up-to-date."
