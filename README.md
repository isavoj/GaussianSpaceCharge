# GaussianSpaceCharge

This repository contains the implementation, with some bugs,  of an Analytic Space-Charge Model for Gaussian Beams with cross-plane coupling, based on the research paper by M. Holz and V. Ziemann ([link to the paper](https://uu.diva-portal.org/smash/get/diva2:1160961/FULLTEXT01.pdf)).

## Installation

To install the package, follow the instructions:

```bash
$ git clone https://github.com/isavoj/GaussianSpaceCharge.git
$ cd GaussianSpaceCharge
```

(Optional step, but recommended:) 
Create and activate a virtual environment:
```bash
$ python -m venv env
$ source env/bin/activate  # On Windows, use 'env\Scripts\activate'
```
Now, 
```bash
$ python setup.py install
```

## Run the main script:
```bash
$ cd GaussianSpaceCharge
$ (env) python main.py
```

 ----------------------------------------
## Project Overview

The project contains 1 module called `GaussianSpaceCharge`, within which you'll find 4 files:
- `main.py`
- `space_charge_calc.py`
- `elements.py`
- `beam.py`

When you run `main.py`, the `main()` function is invoked, which in turn calls the function `propagate_beam_through_lattice()`. This function iterates over each sliced element (divided by 2 to apply Space Charge (SC) more accurately) in a series of FODO cells and calculates the space charge effects for each slice.

### Main Steps in `main()` Function

The process within the `main()` function includes the following key steps:
- **Defining the Lattice Configuration and Properties**: Setup the layout and parameters of the lattice elements used in the simulation.
- **Setting Up the Twiss Parameters for the Beam**: (alpa_x, beta_x, eps_x + y )
- **Configuring the Standard Deviation (SD) and Cutoff for the Gaussian Distribution**: These parameters are used to initialize phase space coordinates.
- **Calculating the Initial Covariance Matrix (`sigma`)**: Derived from the phase space coordinates to represent the beam's initial state.
- **Propagating the Beam Through the Lattice**: The beam is moved through the lattice, applying space charge effects at each step.
- **Plotting the Resulting Beam Envelope**: Visual representation of the beam's evolution throughout the lattice.

### Customization Options

To customize other parts of the simulation:
1. **Beam Parameters**: Change beam parameters like perveance by modifying the `beam_perveance` function in `beam.py`.
2. **Distribution Adjustments**: Adjust the Gaussian and KV distributions directly within `beam.py`.

----------------------------------------



