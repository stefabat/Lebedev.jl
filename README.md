# Lebedev.jl

[![Build Status](https://github.com/stefabat/Lebedev.jl/workflows/CI/badge.svg)](https://github.com/stefabat/Lebedev.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/stefabat/Lebedev.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/stefabat/Lebedev.jl)

A julia package to compute the [Lebedev quadrature rule](https://www.wikiwand.com/en/Lebedev_quadrature)
over the surface of a three-dimensional unit sphere.

This formula allows to approximate a surface integral as follows

<img src="https://latex.codecogs.com/svg.latex?I(f)=\int_0^{2\pi}d\phi\int_0^\pi&space;f(\phi,\theta)\sin\theta&space;d\theta\approx4\pi\sum_iw_if(x_i,y_i,z_i)"/>

where the number of points in the sum depends on the quadrature order chosen.
The coefficients generated by the Lebedev rule allow to calculate integrals of polynomials

<img src="https://latex.codecogs.com/svg.latex?x^ky^lz^m,\;k&plus;l&plus;m\le131"/>

with a relative accuracy in the order of 5⋅10⁻¹⁴. The points created by the rule
have octahedral rotation (Oₕ point group) and inversion symmetry, while the coefficients are
normalized such that they sum up to one.

## Installation

The `Lebedev.jl` package is tested with Julia 1.5 and 1.6 on Linux, MacOS and Windows. Previous versions
of Julia starting at 1.0 might be supported as well, but are not tested; use it at your own risk.

To install the latest version of the package, simply enter in the Julia package manager by typing `]`
in the REPL and issue
```julia
pkg> add Lebedev
```
That's it! You're good to go.

## Usage

After installing the package, type
```julia
julia> using Lebedev
```
in the REPL or add it to your script.

The `Lebedev.jl` package provides essentially two functions to generate quadrature points and weights:
```julia
lebedev_by_points(n::Integer) -> x,y,z,w
lebedev_by_order(n::Integer) -> x,y,z,w
```
Both of them return four one-dimensional arrays of the same length, the first three containing
the `x`, `y` and `z` Cartesian coordinates of the quadrature points (which lie on the unit sphere)
and the last one containing the associated weights `w`.
The `lebedev_by_points` function returns the points and weights corresponding to the `n`-point Lebedev
rule, while the `lebedev_by_order` function returns the points and weights ensuring "exact"
integration for a polynomial of order up to `n`.

Lebedev determined 65 quadrature rules ranging from order 3 up to order 131, increasing two by two
(hence only odd-numbered orders are known), however, this package only provides a subset of them.
Use the function `isavailable(n::Integer)` to know if the quadrature rule for a given order `n` is
available and the function `availablerules()` to print out all available orders and corresponding
number of points. If you need to access the list of available orders or points, use
`getavailablerules()` and `getavailablepoints()`, which return lists sorted in ascending order.

## Examples

```julia
julia> @time lebedev_by_points(2702);
  0.000053 seconds (8 allocations: 84.812 KiB)

julia> @time lebedev_by_order(89);
  0.000035 seconds (8 allocations: 84.812 KiB)

julia> f(x,y,z) = x^2 * y^4 * z^6
f (generic function with 1 method)

# we need at least a rule of order n = 2+4+6 = 12
julia> isavailable(12)
false

# only odd-numbered rules are implemented, check if n = 13 is available
julia> isavailable(13)
true

julia> @time x,y,z,w = lebedev_by_order(13);
  0.000007 seconds (4 allocations: 2.625 KiB)

# integrates f(x,y,z) = x²y⁴z⁶ on the unit sphere
julia> @time 4 * pi * dot(w,f.(x,y,z))
  0.000063 seconds (5 allocations: 768 bytes)
0.0041846055991871965
```
## Reference

Vyacheslav Lebedev, Dmitri Laikov,
“A quadrature formula for the sphere of the 131st algebraic order of accuracy”,
*Doklady Mathematics*, **59** (3), 477-481 (1999).

## License

The Lebedev.jl package is released under the GNU General Public License, version 3.0.
This implementation is based on a C source code developed by Dmitri Laikov and John Burkardt,
which can be found at
https://people.sc.fsu.edu/~jburkardt/c_src/sphere_lebedev_rule/sphere_lebedev_rule.html.
