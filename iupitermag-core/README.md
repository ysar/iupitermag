# iupitermag-core

This is the Rust core of the [`iupitermag`](https://github.com/ysar/iupitermag) Python package, 
also usable as Rust crate.

This library is used to model the magnetic field environment around Jupiter, accounting for both the 
internal field and the current-sheet field. In addition to using named field models like JRM33, 
JRM09, CON2020, you may also use your own model parameters to define a custom field model.

The documentation here is specific to the Rust core, if you'd like to use it directly. Documentation 
on the Python package is available on [Github](https://github.com/ysar/iupitermag). The Rust API
is also documented to some extent on `docs.rs`.

![Jupiter's Surface Field Strength](https://raw.githubusercontent.com/ysar/iupitermag/refs/heads/main/images/jupiter_surfacefield.png)

## Installation
To add this core library in your Rust projects, use
```shell
cargo add iupitermag
```
Or add it to the `Cargo.toml` file.


## Usage 

There are some structs and methods that you can use. Most have an equivalent in
the downstream Python package.

### Calculating the internal and current sheet fields at a single point.

All positions should be in the IAU_JUPITER coordinate system.

```rust
use iupitermag::internal::InternalField;
use iupitermag::currentsheet::{CurrentSheetField, IntegrationType};
use iupitermag::field::Field; // This is the Field trait, necessary for calc_field

// Instantiate the JRM33 and CON2020 fields.
let internal_field = InternalField::new("JRM33", None, None, None);
let currentsheet_field = CurrentSheetField::new(
    "CON2020".to_string(), None, IntegrationType::Analytic
);

let r: f64 = 10.;
let theta: f64 = 0.;
let phi: f64 = 0.;
let b_int_rtp = internal_field.calc_field(r, theta, phi);
let b_ext_rtp = currentsheet_field.calc_field(r, theta, phi);

// Alternatively, you can use the cartesian form
let x: f64 = 10.;
let y: f64 = 0.;
let z: f64 = 0.;
let b_int_xyz = internal_field.calc_field_xyz(x, y, z);
let b_ext_xyz = currentsheet_field.calc_field_xyz(x, y, z);
```

### Calculating the internal and current sheet fields for a collection of points.

If you have a collection of points stored as `ndarray::Array2` of shape (N, 3), 
you can use `map_calc_field` or `parmap_calc_field` (or their corresponding 
cartesian versions `map_calc_field_xyz` and `parmap_calc_field_xyz`).

### Tracing magnetic field lines

`iupitermag` can trace magnetic field lines to Jupiter using `trace_field_to_planet`, which takes 
as input a collection of starting points (which each result in a separate trace).  The coordinates 
for these starting points should be cartesian.

```rust
use iupitermag::currentsheet::{CurrentSheetField, IntegrationType};
use iupitermag::internal::InternalField;
use iupitermag::trace::trace_field_to_planet;
use ndarray::Array;

let internal_field = InternalField::new("JRM33", None, None, Some(10));

let currentsheet_field =
    CurrentSheetField::new("CON2020".to_string(), None, IntegrationType::Analytic);

let start_position = Array::from_vec(vec![-10.0, 2.0, 3.0]);

let trace = trace_field_to_planet(start_position, &internal_field, &currentsheet_field);
// `trace` is Array2 with shape (N, 3) where N is the number of points in the trace and 3 
// refers to the cartesian coordinates of each point.
```

![Traced field lines](https://raw.githubusercontent.com/ysar/iupitermag/refs/heads/main/images/traced_field_lines.png)
