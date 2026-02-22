"""
Methods to calculate field quantities like the JRM magnetic field
model and the CON2020 current sheet model.
"""

import numpy as np


def calc_assoc_legendre_poly(theta, N):
    """
    Calculates Gaussian normalized Legendre polynomials
    Assumes that the input is of the form P_nm(cos(theta)) NOT P_nm(x)
    """
    sintheta = np.sin(theta)
    costheta = np.cos(theta)

    P = np.zeros((N + 1, N + 1))
    K = np.zeros((N + 1, N + 1))
    diffP = np.zeros((N + 1, N + 1))

    P[0, 0] = 2
    diffP[0, 0] = 0

    n = 1
    P[1, 1] = sintheta * P[n - 1, n - 1]
    diffP[1, 1] = sintheta * diffP[n - 1, n - 1] + costheta * P[n - 1, n - 1]
    P[1, 0] = costheta * P[0, 0]
    diffP[1, 0] = costheta * diffP[0, 0] - sintheta * P[0, 0]

    for n in range(2, N + 1):
        P[n, n] = sintheta * P[n - 1, n - 1]
        diffP[n, n] = sintheta * diffP[n - 1, n - 1] + costheta * P[n - 1, n - 1]

        for m in range(0, n):
            K[n, m] = ((n - 1) ** 2 - m**2) / ((2 * n - 1) * (2 * n - 3))
            P[n, m] = costheta * P[n - 1, m] - K[n, m] * P[n - 2, m]
            diffP[n, m] = (
                costheta * diffP[n - 1, m]
                - sintheta * P[n - 1, m]
                - K[n, m] * diffP[n - 2, m]
            )

    return P, diffP


def calc_internal_field(r_a, theta, phi, g, h, N):
    """
    Calculates Psi, Br, Btheta, Bphi based on the spherical harmonics
    representation.
    """
    Br = 0
    Btheta = 0
    Bphi = 0

    a_r = 1.0 / r_a

    P, diffP = calc_assoc_legendre_poly(theta, N)

    sintheta = np.sin(theta)
    sinphi = np.sin(phi)
    cosphi = np.cos(phi)
    inv_sintheta = 0.0
    if np.abs(sintheta) > 1.0e-3:
        inv_sintheta = 1.0 / sintheta

    for n in range(0, N + 1):
        m = 0
        sinphi_prev = 0.0
        cosphi_prev = 1.0
        # Psi += a * (a_r)**(n + 1) * P[n, m] * g[n, m]
        Br += (a_r) ** (n + 2) * (n + 1) * P[n, m] * g[n, m]
        Btheta -= (a_r) ** (n + 2) * diffP[n, m] * g[n, m]
        Bphi += inv_sintheta * (a_r) ** (n + 2) * P[n, m] * m * (-h[n, m])

        for m in range(1, n + 1):
            sinmphi = sinphi_prev * cosphi + cosphi_prev * sinphi
            cosmphi = cosphi_prev * cosphi - sinphi_prev * sinphi

            Br += (
                (a_r) ** (n + 2)
                * (n + 1)
                * P[n, m]
                * (g[n, m] * cosmphi + h[n, m] * sinmphi)
            )

            Btheta -= (
                (a_r) ** (n + 2) * diffP[n, m] * (g[n, m] * cosmphi + h[n, m] * sinmphi)
            )

            Bphi += (
                inv_sintheta
                * (a_r) ** (n + 2)
                * P[n, m]
                * m
                * (g[n, m] * sinmphi - h[n, m] * cosmphi)
            )

            sinphi_prev = sinmphi
            cosphi_prev = cosmphi

        B = np.array([Br, Btheta, Bphi])

    return B


def read_coefficients_file(filename, degree=None, doNormalize=False):
    """Reads the coefficients and return after normalization."""

    N = degree if degree else 30
    g = np.zeros((N + 1, N + 1))
    h = np.zeros((N + 1, N + 1))

    with open(filename, "r") as f:
        f.readline()
        for n in range(1, N + 1):
            for m in range(0, n + 1):
                g[n, m] = float(f.readline().split()[1])
                # print('g', n, m, g[n,m])

            for m in range(1, n + 1):
                h[n, m] = float(f.readline().split()[1])
                # print('h', n, m, h[n,m])

    # Normalize the Gauss coefficients
    if doNormalize:
        S = np.zeros((N + 1, N + 1))
        S[0, 0] = 1
        for n in range(1, N + 1):
            S[n, 0] = S[n - 1, 0] * (2 * n - 1) / n
            S[n, 1] = S[n, 0] * np.sqrt(n * 2.0 / (n + 1))

            g[n, 0] = g[n, 0] * S[n, 0]
            h[n, 0] = h[n, 0] * S[n, 0]
            g[n, 1] = g[n, 1] * S[n, 1]
            h[n, 1] = h[n, 1] * S[n, 1]

            for m in range(2, n + 1):
                S[n, m] = S[n, m - 1] * np.sqrt((n - m + 1) / (n + m))
                g[n, m] = g[n, m] * S[n, m]
                h[n, m] = h[n, m] * S[n, m]

    return g, h, N
