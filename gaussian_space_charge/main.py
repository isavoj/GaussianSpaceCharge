import numpy as np
import cmath
from matplotlib import pyplot as plt
from beam import Beam
from elements import Elements
import space_charge_calc

""" INSTRUCTIONS

When you run this script, the `main()` function is called, which in turn invokes the important function 
`propagate_beam_through_lattice()`.
 This function iterates over each sliced element (quadrupole and drift) 
 and calculates the space charge for each slice.


1. To change beam parameters, modify the `beam_perveance` function in `beam.py`.
2. The creation of Gaussian and KV distributions is also in the `beam.py` script.
3. Twiss parameters are adjusted in the main function of this script.
4. To adjust the standard deviation and cutoff value of the Gaussian distribution, 
modify the `factor` and `cut_off` parameters in the main function of this script.

"""

def calculate_covariance_matrix(x_beam_coordinates, y_beam_coordinates):
    # Extract individual coordinate arrays
    x_values, x_p_values = zip(*x_beam_coordinates)
    y_values, y_p_values = zip(*y_beam_coordinates)
    all_coordinates = np.vstack((x_values, x_p_values, y_values, y_p_values))

    return np.cov(all_coordinates)

def propagate_beam_through_lattice(sigma, lattice_dict, lattice, split_quad, L_quad, split_drift, L_drift):
    sigmax = []
    sigmay = []
    positions = []
    inc_len = 0

    #Go through each element
    for ele in lattice:

        if ele == "D":
            ds = L_drift / (2 * split_drift)
            loop_range = split_drift
        else:
            ds = L_quad / (2 * split_quad)
            loop_range = split_quad

        # Go through each slice of the element
        for _ in range(loop_range):
            inc_len += ds * 2
            positions.append(inc_len)

            # Get external field transfer matrix
            R = lattice_dict[ele]

            # Get SC - inhomogenous terms
            T_l = space_charge_calc.calculate_matrix_T(sigma, ds * 2)

            # Calculate transport half way W.O SC
            sigma = np.dot(np.dot(R, sigma), R.T)

            # Apply SC
            sigma = sigma + T_l

            # Calculate transport half way W SC
            sigma = np.dot(np.dot(R, sigma), R.T)

            sigmax.append(cmath.sqrt(sigma[0, 0]) * 1e3) #x position
            sigmay.append(cmath.sqrt(sigma[2, 2]) * 1e3) #y position

    return positions, sigmax, sigmay

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

    # Parameters
    alpha_x = 0.0
    beta_x = 19.8176
    eps_x = 1.25 * 1e-6

    alpha_y = 0.0
    beta_y = 5.672
    eps_y = 2.5 * 1e-6

    num_particles = int(1e5) # Macro particles

    factor = 1 # This is essentially sigma for the gaussian distribution / (np.sqrt(2)) #(2 * 2.3548) #1 / np.sqrt(2)
    cut_off = 3.1
    split_quad = 20 #nbr of splits
    L_quad = 0.5
    split_drift = 20
    L_drift = 5

    # Initialize coordinates Gaussian
    x_beam_coordinates = Beam.get_coordinates_gauss(alpha_x, beta_x, eps_x, num_particles, factor, cut_off)
    y_beam_coordinates = Beam.get_coordinates_gauss(alpha_y, beta_y, eps_y, num_particles, factor, cut_off)

    # Initialize coordinates KV for comparison
    #x_beam_coordinates, y_beam_coordinates = Beam.get_coordinates_kv(alpha_x, beta_x, eps_x, alpha_y, beta_y, eps_y, num_particles)

    # Calculate covariance matrix
    sigma = calculate_covariance_matrix(x_beam_coordinates, y_beam_coordinates)

    # Define lattice elements
    QD = Elements.quadrupole(L_quad / (2 * split_quad), -0.2)
    QF = Elements.quadrupole(L_quad / (2 * split_quad), 0.2)
    D = Elements.drift(L_drift / (2 * split_drift))

    lattice_dict = {"D": D, "QF": QF, "QD": QD}
    lattice = ["QF", "D", "QD", "QD", "D", "QF"] * 4

    # Propagate beam through lattice
    positions, sigmax, sigmay = propagate_beam_through_lattice(sigma, lattice_dict, lattice, split_quad, L_quad, split_drift, L_drift)

    # Plot results
    plot_results(positions, sigmax, sigmay, factor)

if __name__ == "__main__":
    main()
