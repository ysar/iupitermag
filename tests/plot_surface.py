import numpy as np
import iupitermag as im
import JupiterMag as jm
from matplotlib import pyplot as plt

def calc_field_jupitermag(pos):

    posOut = np.zeros(pos.shape)
    jm.Internal.Config(Model="jrm09",CartesianIn=True,CartesianOut=True)
    for i in range(pos.shape[0]):
        _val = jm.Internal.Field(pos[i, 0], pos[i, 1], pos[i, 2])
        posOut[i, 0], posOut[i, 1], posOut[i, 2] = _val

    return posOut


def calc_field_iupitermag(pos):
    internalField = im.InternalField("JRM09")
    return internalField.map_calc_field(pos)

# -----------------------------------------------------------------------------
f = 1/15.4
a = 1.
b = a
c = (1 - f) * a 

theta_arr = np.linspace(0, np.pi, 100)
phi_arr = np.linspace(0, 2 * np.pi, 100)

theta, phi = np.meshgrid(theta_arr, phi_arr, indexing='ij')

theta = np.ravel(theta)
phi = np.ravel(phi)

posxyz = np.zeros((100 * 100, 3))
posxyz[:, 0] = a * np.sin(theta) * np.cos(phi)
posxyz[:, 1] = b * np.sin(theta) * np.sin(phi)
posxyz[:, 2] = c * np.cos(theta)

pos = np.zeros((100 * 100, 3))
pos[:, 0] = np.sqrt(np.sum(posxyz**2, axis=1))
pos[:, 1] = np.arccos(posxyz[:, 2] / pos[:, 0])
pos[:, 2] = np.arctan2(posxyz[:, 1], posxyz[:, 0])

magArray = calc_field_iupitermag(pos)
# magArray[np.isnan(magArray)] = 0.
for row in magArray[:5]:
    print(row)
Bmag = np.sqrt(np.sum(magArray**2, axis=1)).reshape((100, 100))

fig = plt.figure(dpi=300, figsize=(4, 3))
ax = plt.subplot(111)
mp=ax.contourf(phi_arr * 180 / np.pi, theta_arr * 180 / np.pi, Bmag, cmap='inferno')
ax.set_title(r'$|\mathbf{B}|$') # 'Using iupiterMag')
ax.set_xlabel(r"$\phi$ [$^\circ$]")
ax.set_ylabel(r"$\theta$ [$^\circ$]")
cbar = fig.colorbar(mp)
cbar.set_label("nT")
fig.savefig("jupiter_surfacefield.png", facecolor='w', bbox_inches='tight')

plt.show()
plt.close(fig)


# magArray1 = calc_field_iupitermag(pos)
# Bmag1 = np.sqrt(np.sum(magArray1**2, axis=1)).reshape((100, 100))
# magArray2 = calc_field_jupitermag(pos)
# Bmag2 = np.sqrt(np.sum(magArray2**2, axis=1)).reshape((100, 100))

# fig = plt.figure(dpi=300)
# ax = plt.subplot(111)
# mp=ax.contourf(phi_arr * 180 / np.pi, theta_arr * 180 / np.pi, Bmag2 - Bmag1, cmap='bwr')
# maxR = np.abs((Bmag2 - Bmag1) / Bmag2).max()
# ax.set_title(f'Difference between other models\n' + \
#              f'Max relative |error| = {maxR * 100:.2e} %\n' + \
#              f"Max absolute |error| = {maxR:.2e} nT"   )
# ax.set_xlabel(r"$\phi$ [deg]")
# ax.set_ylabel(r"$\theta$ [deg]")
# cbar = fig.colorbar(mp)
# cbar.set_label("nT")

# fig.savefig("jupiter_surfacefield_diff.png", facecolor='w', bbox_inches='tight')

plt.show()
plt.close(fig)