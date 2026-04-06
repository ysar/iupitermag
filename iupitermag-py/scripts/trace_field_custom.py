import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

import iupitermag as im

plt.style.use("dark_background")


def plot_jupiter(ax, **kwargs):
    jupiter = Ellipse((0.0, 0.0), 2.0, 2.0 * (1.0 - 1.0 / 15.4), **kwargs)
    ax.add_patch(jupiter)


if __name__ == "__main__":
    # Create a new axially aligned dipole field using Schmidt coefficients.
    internal_field = im.InternalField(
        "Custom",
        g=np.array([[0.0, 0.0], [410993.4, 0.0]]),
        h=np.array([[0.0, 0.0], [0.0, 0.0]]),
    )

    # Get parameters for CON2020
    con2020_field = im.CurrentSheetField("CON2020")
    params = con2020_field.get_params()

    # Change current sheet tilt to 0.0
    params["theta_d"] = 0.0

    # Define a modified current sheet field based on CON2020 parameters.
    currentsheet_field = im.CurrentSheetField("Custom", params=params)

    starting_positions = np.array(
        [
            [-10.0, 0.0, 0.0],
            [-15.0, 0.0, 0.0],
            [-20.0, 0.0, 0.0],
            [-25.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [15.0, 0.0, 0.0],
            [20.0, 0.0, 0.0],
            [25.0, 0.0, 0.0],
        ]
    )

    trace = im.trace_field_to_planet(
        starting_positions, internal_field, currentsheet_field
    )

    fig, ax = plt.subplots(1, 1, figsize=(6, 3), dpi=200)
    for t in trace:
        ax.plot(t[:, 0], t[:, 2], lw=1.0)
    ax.axhline(0.0, color="0.5", lw=0.25)
    ax.axvline(0.0, color="0.5", lw=0.25)
    plot_jupiter(ax, color="0.8", zorder=3)
    ax.set_xlabel(r"X [R$_J$]")
    ax.set_ylabel(r"Z [R$_J$]")
    ax.set_aspect("equal")
    fig.savefig(
        "../images/traced_field_lines_custom.png", facecolor="k", bbox_inches="tight"
    )
    plt.show()
    plt.close(fig)
