# Internal Field

The internal magnetic field of the planet body can be described via spherical harmonics. E.g. the scalar magnetic potential \( A \) is,
 
\[ A(r, \theta, \phi) = a \sum_{n=1}^k \left(\frac{a}{r}\right)^{l+1} \sum_{m=0}^{n} [g_n^m \cos(m\phi) + h_n^m \sin(m\phi)] P_n^m (\cos\theta)  \]

From \(A \), we can calculate the components of magnetic field intensity (\( B_r, B_\theta, B_\phi \))
at any point. \(g_n^m\) and \(h_n^m\) are coefficients (typically in nT) defined
for some degree \(N \), where \(n = 1, 2, 3, ..., N\) and \(m = 0, 1, 2, ... n \).

Before using an internal field model, you need to instantiate an `InternalField`
class object, like so.

```python
import iupitermag as iu

internalField = iu.InternalField("JRM09")
```

You can choose to truncate a field to a specified degree. 
```python
internalField = iu.InternalField("JRM09", degree=5)
```

You can also create your own field by passing in arrays representing the coefficients.
`iupitermag` defines these arrays such that `g[1, 2]` corresponds to \(g_1^2\). So,
arrays will have a zeroth row corresponding to \( n=0 \), that should be set to zero.
 Coefficients arrays should have shape `(degree + 1, degree + 1)`.
```python
# This will create an axially aligned dipole field
g = np.array([[0., 0.], [21000., 0.]])
h = np.array([[0., 0.], [0., 0.]])
internalField = iu.InternalField("Custom", g=g, h=h, degree=1)
```

Once you have defined a field, you can call various methods on it. 

**Calculate the magnetic field vector at a specific point using `calc_field`**
```python
r = 10.
theta = 0.5 * np.pi
phi = 0.
b = internalField.calc_field(r, theta, phi)
```

**Calculate the magnetic field vector at many points using `map_calc_field`**
```python
points = np.zeros((1000, 3))
points[:, 0] = 10.    # Set the radius to 10
bPoints = internalField.map_calc_field(points)
```