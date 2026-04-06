import numpy as np

import iupitermag as im


def calc_field_iupitermag(pos):
    internalField = im.InternalField("JRM33")
    return internalField.map_calc_field(pos)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    f = 1 / 15.4
    a = 1.0
    b = a
    c = (1 - f) * a

    N = 1000
    theta_arr = np.linspace(0, np.pi, N)
    phi_arr = np.linspace(0, 2 * np.pi, N)

    theta, phi = np.meshgrid(theta_arr, phi_arr, indexing="ij")

    theta = np.ravel(theta)
    phi = np.ravel(phi)

    pos_xyz = np.zeros((N * N, 3))
    pos_xyz[:, 0] = a * np.sin(theta) * np.cos(phi)
    pos_xyz[:, 1] = b * np.sin(theta) * np.sin(phi)
    pos_xyz[:, 2] = c * np.cos(theta)

    pos = np.zeros((N * N, 3))
    pos[:, 0] = np.sqrt(np.sum(pos_xyz**2, axis=1))
    pos[:, 1] = np.arccos(pos_xyz[:, 2] / pos[:, 0])
    pos[:, 2] = np.arctan2(pos_xyz[:, 1], pos_xyz[:, 0])

    mag_array = calc_field_iupitermag(pos)
