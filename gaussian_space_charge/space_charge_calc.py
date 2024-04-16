import numpy as np
import cmath
import math
from beam import Beam


"""Class for NonLinear space charge calculations and 
two standalone functions for Linear space charge 
and calculation of the T matrix."""
class NonLinearSc:

    @staticmethod
    def S11(sigma, n, c):
        term1 = 1.0
        term2 = 2.0 / n * (sigma[0, 0] * c[0] ** 2 + sigma[2, 2] * c[1] ** 2 + 2 * sigma[0, 2] * c[0] * c[1])
        return term1 + term2

    @staticmethod
    def S12(sigma, n, c, d):
        term1 = 2.0 / n * sigma[0, 0] * c[0] * d[0] + sigma[2, 2] * c[1] * d[1] + sigma[0, 2] * (c[0] * d[1] + c[1] * d[0])
        return term1

    @staticmethod
    def S22(sigma, n, d):
        term1 = 1.0
        term2 = 2.0 / n * (sigma[0, 0] * d[0] ** 2 + sigma[2, 2] * d[1] ** 2 + 2 * sigma[0, 2] * d[0] * d[1])
        return term1 + term2

    @staticmethod
    def z1(sigma):
        denom = cmath.sqrt(2 * (sigma[0, 0] - sigma[2, 2] + 2j * sigma[0, 2]))
        a1 = 1 / denom
        a3 = 1j / denom
        denom_bar = cmath.sqrt(2 * (sigma[0, 0] - sigma[2, 2] - 2j * sigma[0, 2]))
        a1_bar = 1 / denom_bar
        a3_bar = -1j / denom_bar
        return a1, a3, a1_bar, a3_bar

    @staticmethod
    def z2(sigma):
        term1 = cmath.sqrt(sigma[0, 0] * sigma[2, 2] - sigma[0, 2] ** 2)
        term2 = cmath.sqrt(2 * (sigma[0, 0] - sigma[2, 2] + 2j * sigma[0, 2]))
        denom = term1 * term2
        b1 = (sigma[2, 2] - 1j * sigma[0, 2]) / denom
        b3 = (1j * sigma[0, 0] - sigma[0, 2]) / denom
        denom_bar = term1 * cmath.sqrt(2 * (sigma[0, 0] - sigma[2, 2] - 2j * sigma[0, 2]))
        b1_bar = (sigma[2, 2] + 1j * sigma[0, 2]) / denom_bar
        b3_bar = (-1j * sigma[0, 0] - sigma[0, 2]) / denom_bar

        return b1, b3, b1_bar, b3_bar

    @classmethod
    def P(cls,sigma, n, c, d):
        S_11 = cls.S11(sigma, n, c)
        S_12 = cls.S12(sigma, n, c, d)
        S_22 = cls.S22(sigma, n, d)

        inner_term = np.pi / 2 - np.arctan(S_12 / cmath.sqrt(S_11 * S_22 - S_12 ** 2))
        return 2 / (n * np.pi * cmath.sqrt(S_11 * S_22 - S_12 ** 2)) * inner_term

    @classmethod
    def calc_F0F0(cls,sigma):
        a1, a3, a1_bar, a3_bar = NonLinearSc.z1(sigma)
        b1, b3, b1_bar, b3_bar = NonLinearSc.z2(sigma)
        P1 = cls.P(sigma, 1, [a1, a3], [a1_bar, a3_bar])
        P2 = cls.P(sigma, 2, [a1, a3], [b1_bar, b3_bar])
        P3 = cls.P(sigma, 2, [a1_bar, a3_bar], [b1, b3])
        P4 = cls.P(sigma, 3, [b1, b3], [b1_bar, b3_bar])
        term1 = -math.pi / cmath.sqrt(2 * (sigma[0, 0] - sigma[2, 2] + 2j * sigma[0, 2]))
        term2 = 1 / cmath.sqrt(2 * (sigma[0, 0] - sigma[2, 2] - 2j * sigma[0, 2]))
        FOFO = term1 * (P1 - P2 - P3 + P4) * term2
        return FOFO


def calculate_T_components_linear( sigma_k1, sigma_k3, sigma11, sigma33, sigma13):
    denominator = 2 * (sigma11 - sigma33 + 2j * sigma13)
    inner_term_numerator = sigma_k1 * (sigma33 - 1j * sigma13) + 1j * sigma_k3 * (sigma11 + 1j * sigma13)
    inner_term_denominator = cmath.sqrt(sigma11 * sigma33 - sigma13 ** 2)
    result = (1j / denominator) * (sigma_k1 + 1j * sigma_k3 - inner_term_numerator / inner_term_denominator)
    return result.imag, result.real

def calculate_matrix_T(covariance_matrix, ds):

        sigma11 = covariance_matrix[0, 0]
        sigma33 = covariance_matrix[2, 2]
        sigma13 = covariance_matrix[0, 2]
        sigma21 = covariance_matrix[1, 0]
        sigma23 = covariance_matrix[1, 2]
        sigma31 = covariance_matrix[2, 0]
        sigma41 = covariance_matrix[3, 0]
        sigma43 = covariance_matrix[3, 2]

        K = Beam.beam_perveance() * ds

        # Linear components
        x1_f1, x1_f3 = calculate_T_components_linear(sigma11, sigma13, sigma11, sigma33, sigma13)  # Real
        x2_f1, x2_f3 = calculate_T_components_linear(sigma21, sigma23, sigma11, sigma33, sigma13)  # Real
        x3_f1, x3_f3 = calculate_T_components_linear(sigma31, sigma33, sigma11, sigma33, sigma13)
        x4_f1, x4_f3 = calculate_T_components_linear(sigma41, sigma43, sigma11, sigma33, sigma13)

        # Quadratic components
        FOFO = NonLinearSc.calc_F0F0(covariance_matrix)
        f3_squared = FOFO / 2
        # f1_squared = FOFO.imag / 2
        f1_squared = FOFO / 2
        f1_f3 = 0

        # WITH the x2_f1 term and the x4_f3 term.

        # T = np.array([
        #     [0, K * x1_f1, 0, K * x1_f3],
        #     [0, K * 2 * x2_f1 + K ** 2 * f1_squared, K * x3_f1, K ** 2 * f1_f3 + K * x2_f3 + K * x4_f1],
        #     [0, 0, 0, K * x3_f3],
        #     [0, 0, 0, K * 2 * x4_f3 + K ** 2 * f3_squared]
        # ])


        # WITHOUT the x2_f1 term and the x4_f3 term.
        T = np.array([
            [0, K * x1_f1, 0, K * x1_f3],
            [0, K ** 2 * f1_squared, K * x3_f1, K ** 2 * f1_f3 + K * x2_f3 + K * x4_f1],
            [0, 0, 0, K * x3_f3],
            [0, 0, 0, K ** 2 * f3_squared]
        ])

        return T
