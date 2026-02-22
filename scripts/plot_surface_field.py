import numpy as np
from matplotlib import pyplot as plt

import iupitermag as im

plt.style.use("dark_background")


def calc_field_iupitermag(pos):
    internalField = im.InternalField("JRM33")
    return internalField.map_calc_field(pos)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    f = 1 / 15.4
    a = 1.0
    b = a
    c = (1 - f) * a

    theta_arr = np.linspace(0, np.pi, 100)
    phi_arr = np.linspace(0, 2 * np.pi, 100)

    theta, phi = np.meshgrid(theta_arr, phi_arr, indexing="ij")

    theta = np.ravel(theta)
    phi = np.ravel(phi)

    pos_xyz = np.zeros((100 * 100, 3))
    pos_xyz[:, 0] = a * np.sin(theta) * np.cos(phi)
    pos_xyz[:, 1] = b * np.sin(theta) * np.sin(phi)
    pos_xyz[:, 2] = c * np.cos(theta)

    pos = np.zeros((100 * 100, 3))
    pos[:, 0] = np.sqrt(np.sum(pos_xyz**2, axis=1))
    pos[:, 1] = np.arccos(pos_xyz[:, 2] / pos[:, 0])
    pos[:, 2] = np.arctan2(pos_xyz[:, 1], pos_xyz[:, 0])

    mag_array = calc_field_iupitermag(pos)

    Bmag = np.sqrt(np.sum(mag_array**2, axis=1)).reshape((100, 100))

    fig = plt.figure(dpi=200, figsize=(6, 3))
    ax = plt.subplot(111, projection="mollweide")

    longitude = phi_arr - np.pi
    latitude = np.pi / 2 - theta_arr

    mp = ax.contourf(longitude, latitude, Bmag, cmap="inferno")
    ax.set_title(r"$|\mathbf{B}|$")
    # ax.set_xlabel(r"Lon [$^\circ$]")
    ax.set_ylabel(r"Lat [$^\circ$]")
    ax.grid(True, color="0.5", lw=0.25)
    ax.set_xticklabels([])

    cbar = fig.colorbar(mp)
    cbar.set_label("nT")
    fig.savefig("images/jupiter_surfacefield.png", facecolor="k", bbox_inches="tight")

    plt.show()
    plt.close(fig)
