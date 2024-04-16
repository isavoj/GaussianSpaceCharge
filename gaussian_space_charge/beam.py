import random
import math
import numpy as np

" Gives you gaussian or kv coordinates, they are normalized but converted into the real-phase space coordinates through the function def get_u_up"
class Beam:

    N = None  # Defining N as a class variable

    @classmethod
    def get_coordinates_kv(cls, alpha_x, beta_x, eps_x, alpha_y, beta_y, eps_y, num_particles):
        """Return coordinates for the 2D KV-distribution."""
        eps_x *= 4
        eps_y *= 4
        coordinates_x = []
        coordinates_y = []

        for _ in range(num_particles):
            phi = 2 * math.pi * (random.random() - 0.5)
            rho = math.sqrt(random.random())
            x_norm = rho * math.cos(phi)
            y_norm = rho * math.sin(phi)

            p0 = math.sqrt(abs(1.0 - rho ** 2))
            phi = 2 * math.pi * (random.random() - 0.5)
            xp_norm = p0 * math.cos(phi)
            yp_norm = p0 * math.sin(phi)

            u_x, up_x, _ = cls.get_u_up(x_norm, xp_norm, alpha_x, beta_x, eps_x)
            u_y, up_y, _ = cls.get_u_up(y_norm, yp_norm, alpha_y, beta_y, eps_y)

            coordinates_x.append([u_x, up_x])
            coordinates_y.append([u_y, up_y])

        return coordinates_x, coordinates_y

    @classmethod
    def get_coordinates_gauss(cls, alpha, beta, emittance, num_particles, factor, cut_off):
        """
        Generate 1D Gaussian distributed coordinates for particles.
        """
        coordinates = []
        sigma = factor
        cut_off2 = cut_off ** 2

        for _ in range(num_particles):
            x_norm = random.gauss(0.0, sigma)
            xp_norm = random.gauss(0.0, sigma)
            if cut_off > 0.0:
                while (x_norm ** 2 + xp_norm ** 2) > cut_off2:
                    x_norm = random.gauss(0.0, sigma)
                    xp_norm = random.gauss(0.0, sigma)

            u, up, _ = cls.get_u_up(x_norm, xp_norm, alpha, beta, emittance)
            coordinates.append([u, up])

        return coordinates

    @classmethod
    def beam_perveance(cls):
        """
        Calculate the perveance of the beam according to paper
        """
        cls.N = 5e13
        e = 1.602176634 * 1e-19
        c = 3e8
        e_kin_ini = 5.0
        m0 = 0.938
        gamma = (m0 + e_kin_ini) / m0
        beta = np.sqrt(gamma * gamma - 1.0) / gamma

        eps0 = 8.854 * 10 ** (-12)
        bunch_size = 0.59205
        cp_massa = 1.6726219e-27
        r0 = e ** 2 / (4 * np.pi * eps0 * cp_massa * c ** 2)

        nominator = cls.N * 2 * r0
        denominator = np.sqrt(2 * np.pi) * bunch_size * beta ** 2 * gamma ** 3
        our_K = nominator / denominator

        return our_K

    @staticmethod
    def get_u_up(u_norm, up_norm, alpha, beta, emittance):
        """
        Convert normalized coordinates to real phase-space coordinates using Twiss parameters.
        """
        u = np.sqrt(beta * emittance) * u_norm
        up = np.sqrt(emittance / beta) * (up_norm - alpha * u_norm)
        return [u, up, 0]
