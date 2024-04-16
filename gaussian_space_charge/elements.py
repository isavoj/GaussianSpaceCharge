import numpy as np

class Elements:

    @staticmethod
    def drift(L: float) -> np.ndarray:
        """
        Creates a drift matrix for a given length 'L'.

        Parameters:
        L (float): Length of the drift.

        Returns:
        numpy.ndarray: 4x4 drift matrix.
        """
        out = np.eye(4)
        out[0, 1] = L
        out[2, 3] = L
        return out

    @staticmethod
    def quadrupole(L: float, k: float) -> np.ndarray:
        """
        Creates a quadrupole matrix for a given length 'L' and quadrupole strength 'k'.

        Parameters:
        L (float): Length of the quadrupole (in meters).
        k (float): Quadrupole strength.

        Returns:
        numpy.ndarray: 4x4 quadrupole matrix.
        """
        D = 1  # np.sqrt(1 + 2 * 0.1 / self.beam.beta + 0.1 ** 2)  # Chromaticity, ignore
        ksq = np.sqrt(np.abs(k) / D)
        out = np.eye(4)

        A = np.array([[np.cos(ksq * L), np.sin(ksq * L) / (ksq * D)],
                      [-ksq * D * np.sin(ksq * L), np.cos(ksq * L)]])
        B = np.array([[np.cosh(ksq * L), np.sinh(ksq * L) / (ksq * D)],
                      [ksq * D * np.sinh(ksq * L), np.cosh(ksq * L)]])

        if k > 0:
            out[0:2, 0:2] = A
            out[2:4, 2:4] = B
        else:
            out[0:2, 0:2] = B
            out[2:4, 2:4] = A

        return out
