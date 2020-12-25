# Lebedev.jl

[![Build Status](https://github.com/stefabat/Lebedev.jl/workflows/CI/badge.svg)](https://github.com/stefabat/Lebedev.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/stefabat/Lebedev.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/stefabat/Lebedev.jl)

A julia package to compute the Lebedev quadrature rule over the surface of a unit sphere.

## References

[1] Axel Becke,
    “A multicenter numerical integration scheme for polyatomic molecules,”
    *Journal of Chemical Physics*, **88** (4), 2547-2553 (1988).

[2] Vyacheslav Lebedev, Dmitri Laikov,
    “A quadrature formula for the sphere of the 131st algebraic order of accuracy,”
    *Russian Academy of Sciences Doklady Mathematics*, **59** (3), 477-481 (1999).

[3] Vyacheslav Lebedev,
    “A quadrature formula for the sphere of 59th algebraic order of accuracy,”
    *Russian Academy of Sciences Doklady Mathematics*, **50**, 283-286 (1995).

[4] Vyacheslav Lebedev, A.L. Skorokhodov,
    “Quadrature formulas of orders 41, 47, and 53 for the sphere,”
    *Russian Academy of Sciences Doklady Mathematics*, **45**, 587-592 (1992).

[5] Vyacheslav Lebedev,
    “Spherical quadrature formulas exact to orders 25-29,”
    *Siberian Mathematical Journal*, **18**, 99-107 (1977).

[6] Vyacheslav Lebedev,
    “Quadratures on a sphere,”
    *Computational Mathematics and Mathematical Physics*, **16**, 10-24 (1976).

[7] Vyacheslav Lebedev,
    “Values of the nodes and weights of ninth to seventeenth order Gauss-Markov quadrature formulae invariant under the octahedron group with inversion,”
    *Computational Mathematics and Mathematical Physics*, **15**, 44-51 (1975).

## License

The Lebedev.jl package is released under the GNU General Public License, version 3.0.
This implementation is based on the C source code originally developed by Dmitri Laikov
and later modified by J. Burkhadt, which can be found at
https://people.sc.fsu.edu/~jburkardt/c_src/sphere_lebedev_rule/sphere_lebedev_rule.html.
