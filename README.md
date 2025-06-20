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

The following example shows how to simulate the 2D Taylor-Green vortex:

1. Clone or download this repository to your computer.
2. Navigate to `sph_code/2d-fixed-wall/source_code/` directory.
3. Run `make` to build the program.
4. Run `./start_calculation` to start the simulation.
5. After the simulation finishes, run `TG_main.py` to generate figures.


The below example is for simulating 2D Taylor-Green vortex
1. Download this code to your PC
2. Move to "./sph_code/2d-fixed-wall/source_code/"
3. Execute `make`
4. Execute `./start_calculation` to simulate
5. After the simulation, execute `TG_main.py` to make figures

## User Manual (Japanese Only)

See the [LS-SPH Fluid Simulator User Manual](./manual.pdf) written in Japanese.  
An English version is not currently available. Please consider using translation tools such as AI.

## License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) for details.


## References

Currently under construction