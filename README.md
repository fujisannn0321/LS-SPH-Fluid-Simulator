# LS-SPH-Fluid-Simulator

This is an open-source fluid simulation code based on the Least Squares Smoothed Particle Hydrodynamics (LS-SPH), a high-precision and generalized extension of classical SPH.

This code provides:
1. A variety of standard fluid benchmark problems, including:
   - 2D Taylor-Green vortex
   - 2D Lid-driven cavity flow
   - 2D Boussinesq convection
2. Applications to geophysical and engineering problems (currently under development)

## Features

- Implementation of both LS-SPH and classical SPH models
- Support for 2D fluid simulations
- Written in Fortran with OpenMP and Python

## Requirements

- Unix-like environment (Linux, macOS, or WSL on Windows)
- [Intel Fortran compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.n7d5f5) (using `ifx`)
- Intel MKL (Math Kernel Library) for linear algebra routines
- `make` utility for building Fortran files
- Python with common libraries (e.g., `matplotlib`, `scipy`)
- `ffmpeg` for movie generation

## Usage

To compile the code:

## User Manual (Japanese Only)

See the [LS-SPH Fluid Simulator User Manual](./manual.pdf) written in Japanese.  
An English version is not currently available. Please consider using translation tools such as AI.

## License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) for details.


## References

Currently under construction
