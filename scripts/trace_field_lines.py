import iupitermag as im
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
plt.style.use("dark_background")


def plot_jupiter(ax, **kwargs):
    jupiter = Ellipse((0.0, 0.0), 2.0, 2.0 * (1.0 - 1.0 / 15.4), **kwargs)
    ax.add_patch(jupiter)
        
if __name__ == "__main__":
    
    internal_field = im.InternalField("JRM33")
    cs_field = im.CurrentSheetField("CON2020")
    
    starting_positions = np.array([
        [-10., 0., 0.],
        [-15., 0., 0.],
        [-20., 0., 0.],
        [-25., 0., 0.],
        [10., 0., 0.],
        [15., 0., 0.],
        [20., 0., 0.],
        [25., 0., 0.],
    ])
    
    trace = im.trace_field_to_planet(starting_positions, internal_field, cs_field)



    fig, ax = plt.subplots(1, 1, figsize=(6,3), dpi=200)
    for t in trace:
        ax.plot(t[:, 0], t[:, 2], lw=1.)
    ax.axhline(0., color='0.5', lw=0.25)
    ax.axvline(0., color='0.5', lw=0.25)
    plot_jupiter(ax, color='0.8', zorder=3)
    ax.set_xlabel(r"X [R$_J$]")
    ax.set_ylabel(r"Z [R$_J$]")
    ax.set_aspect("equal")
    fig.savefig("traced_field_lines.png", facecolor="k", bbox_inches="tight")
    plt.show()
    plt.close(fig)