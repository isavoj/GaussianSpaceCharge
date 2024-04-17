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


## Description

The project contains 1 module called `GaussianSpaceCharge`, within which you'll find 4 files:
- `main.py`
- `space_charge_calc.py`
- `elements.py`
- `beam.py`

### Space- charge 
The space charge calculation code is implemented in `space_charge_calc.py`. The main function `calculate_matrix_T(covariance_matrix, ds)` computes the T matrix, which includes two versions: one with and one without the `x2_f1` and `x4_f3` terms.
comment them out depending on. 
### Execution
Run `main.py`. When you execute this script, the `main()` function is called, which in turn invokes the important function `propagate_beam_through_lattice()`. This function iterates over each sliced element (quadrupole and drift) four our (FODO cell)*4 and calculates the space charge for each slice.

### Modifying additional parameters
1. To change beam parameters, modify the `beam_perveance` function in `beam.py`.
2. The creation of Gaussian and KV distributions is also in the `beam.py` script.
3. Twiss parameters and lattice are adjusted in the main function of this script.
4. To adjust the standard deviation and cutoff value of the Gaussian distribution, modify the `factor` and `cut_off` parameters in the main function of this script.


