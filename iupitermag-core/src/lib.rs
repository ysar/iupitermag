#![doc = include_str!("../README.md")]
#![warn(missing_docs)]

/// Contains methods useful for converting positions and vectors between coordinate systems.
pub mod convert;

/// Definitions for Jupiter's current sheet models.
pub mod currentsheet;

/// Common module for all types of fields. Contains the `Field` trait.
pub mod field;

/// Definitions for Jupiter's internal field models.
pub mod internal;

/// Some functions for calculating Legendre polynomials.
pub mod legendre;

/// Methods for tracing magnetic field lines.
pub mod trace;
