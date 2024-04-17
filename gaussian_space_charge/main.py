import numpy as np
from matplotlib import pyplot as plt
from beam import Beam
from elements import Elements
import space_charge_calc

""" INSTRUCTIONS

Elements: Quads and Drift
Lattice: FODO * nbr_of_cells

----------------------------------------
When you run this script, the `main()` function is invoked, which in turn calls the function
`propagate_beam_through_lattice()`. 
This function iterates over each sliced element( divided by 2 to apply SC more accurately) 
in a series of FODO cells and calculates the space charge effects for each slice.


The main steps in `main()` function include:
- Defining the lattice configuration and properties.
- Setting up the Twiss parameters for the beam.
- Configuring the standard deviation (SD) and cutoff for the Gaussian distribution to initialize phase space coordinates.
- Calculating the initial covariance matrix (`sigma`) from these coordinates.
- Propagating the beam through the lattice.
- Plotting the resulting beam envelope.

To customize other part of the simulation:
1. Change beam parameters like perveance by modifying the `beam_perveance` function in `beam.py`.
2. Adjust the Gaussian and KV distributions directly within `beam.py`.

----------------------------------------
"""

def calculate_covariance_matrix(x_beam_coordinates, y_beam_coordinates):
    # Extract individual coordinate arrays
    x_values, x_p_values = zip(*x_beam_coordinates)
    y_values, y_p_values = zip(*y_beam_coordinates)
    all_coordinates = np.vstack((x_values, x_p_values, y_values, y_p_values))

    return np.cov(all_coordinates)

def propagate_beam_through_lattice(sigma, lattice, lattice_dict):

    sigmax, sigmay, positions = [], [], []
    incremental_length = 0

    for element in lattice:
        element_obj, L_element, split_element = lattice_dict[element]
        ds = L_element / (2 * split_element)  # Calculate slice length
        split_element = int(split_element)  # Ensure split_element is an integer

        for _ in range(split_element):

            incremental_length += ds * 2
            positions.append(incremental_length)

            R = element_obj # Get external field transfer matrix
            T = space_charge_calc.calculate_matrix_T(sigma, ds * 2) # Get SC matrix - inhomogenous terms

            sigma = np.dot(np.dot(R, sigma), R.T) # Calculate transport half way W.O SC

            sigma = sigma + T # Apply SC

            sigma = np.dot(np.dot(R, sigma), R.T) # Calculate transport half way W SC


            sigmax.append(np.sqrt(sigma[0, 0]) * 1e3)
            sigmay.append(np.sqrt(sigma[2, 2]) * 1e3)

    return positions, sigmax, sigmay


def create_FODO_lattice(L_quad, split_quad, L_drift, split_drift, k, cells):
    # Define lattice elements
    QD = Elements.quadrupole(L_quad / (2 * split_quad), -k)
    QF = Elements.quadrupole(L_quad / (2 * split_quad), k)
    D = Elements.drift(L_drift / (2 * split_drift))

    # Create Lattice Dictionary
    lattice_dict = {
        "QD": (QD, L_quad, split_quad),
        "QF": (QF, L_quad, split_quad),
        "D":  (D, L_drift, split_drift)
    }

    # Create Lattice Configuration
    lattice = ["QF", "D", "QD", "QD", "D", "QF"] * cells

    return lattice, lattice_dict

def plot_results(positions, sigmax, sigmay, factor):


    fig, ax1 = plt.subplots(ncols=1, figsize=(18,8))
    ax1.set_title('Gauss with inhomogenous SC-kicks, N = ' + f"{Beam.N:.0e}", fontsize=22)
    ax1.set_xlabel('Position [m]', fontsize=22)
    ax1.set_ylabel(r'$\sigma_x, \sigma_y$ [mm]', fontsize=22)
    ax1.plot(positions, np.divide(sigmax, factor), 'k', label=r'$\sigma_x$')
    ax1.plot(positions, np.divide(sigmay, factor), 'r-.', label=r'$\sigma_y$')
    ax1.legend(fontsize=22)
    plt.show()


def main():

    # Twiss Parameters
    alpha_x = 0.0
    beta_x = 19.8176
    eps_x = 1.25 * 1e-6

    alpha_y = 0.0
    beta_y = 5.672
    eps_y = 2.5 * 1e-6

    # Nbr of Macro Particles
    num_particles = int(1e5) # Macro particles

    # Gaussian Dist parameters
    SD = 1 # This is essentially sigma for the gaussian distribution / (np.sqrt(2)) #(2 * 2.3548) #1 / np.sqrt(2)
    cut_off = 3.1

    # FODO Lattice config
    split_quad = 20  # nbr of splits
    L_quad = 0.5
    split_drift = 20
    L_drift = 5
    k = 0.2
    cells =4 #nbr of FODO cells

    lattice, lattice_dict = create_FODO_lattice(L_quad, split_quad, L_drift, split_drift, k, cells)

    # Initialize coordinates Gaussian
    x_beam_coordinates = Beam.get_coordinates_gauss(alpha_x, beta_x, eps_x, num_particles, SD, cut_off)
    y_beam_coordinates = Beam.get_coordinates_gauss(alpha_y, beta_y, eps_y, num_particles, SD, cut_off)

    # Initialize coordinates KV for comparison
    #x_beam_coordinates, y_beam_coordinates = Beam.get_coordinates_kv(alpha_x, beta_x, eps_x, alpha_y, beta_y, eps_y, num_particles)

    # Calculate covariance matrix
    sigma = calculate_covariance_matrix(x_beam_coordinates, y_beam_coordinates)

    # Propagate beam through lattice
    positions, sigmax, sigmay = propagate_beam_through_lattice(sigma, lattice, lattice_dict)

    # Plot results
    plot_results(positions, sigmax, sigmay, SD)

if __name__ == "__main__":
    main()
