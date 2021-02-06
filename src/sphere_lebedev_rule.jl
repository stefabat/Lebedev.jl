
"""
    gen_oh!(code::Integer, a::Real, b::Real, v::Real, n::Integer,
            x::AbstractArray, y::AbstractArray, z::AbstractArray, w::AbstractArray)

Generate points under Oₕ symmetry.
Given a point on a sphere, specified by `a` and `b`, this function generates
all the equivalent points under Oₕ symmetry, making grid points with weight `v`.
The input arrays `x`, `y`, `z` and `w` are modified in-place starting at index `n+1`.
The function returns the number of elements modified in the input arrays.

Depending on the symmetry group defined by `code`, there are from 6 to 48 different
but equivalent points that are generated:

code=1:   (0,0,1) etc                                  ( 6 points)
code=2:   (0,a,a) etc, a=1/sqrt(2)                     (12 points)
code=3:   (a,a,a) etc, a=1/sqrt(3)                     ( 8 points)
code=4:   (a,a,b) etc, b=sqrt(1 - 2a²)                 (24 points)
code=5:   (a,b,0) etc, b=sqrt(1 - a²), a input         (24 points)
code=6:   (a,b,c) etc, c=sqrt(1 - a² - b²), a, b input (48 points)
"""
function gen_oh!(code::Integer, a::Real, b::Real, v::Real, n::Integer,
                 x::AbstractArray, y::AbstractArray, z::AbstractArray, w::AbstractArray)

    if code == 1
        a = 1.0
        x[n + 1] =   a; y[n + 1] = 0.0; z[n + 1] = 0.0; w[n + 1] = v;
        x[n + 2] =  -a; y[n + 2] = 0.0; z[n + 2] = 0.0; w[n + 2] = v;
        x[n + 3] = 0.0; y[n + 3] =   a; z[n + 3] = 0.0; w[n + 3] = v;
        x[n + 4] = 0.0; y[n + 4] =  -a; z[n + 4] = 0.0; w[n + 4] = v;
        x[n + 5] = 0.0; y[n + 5] = 0.0; z[n + 5] =   a; w[n + 5] = v;
        x[n + 6] = 0.0; y[n + 6] = 0.0; z[n + 6] =  -a; w[n + 6] = v;
        n = 6
    elseif code == 2
        a = 1.0 / sqrt(2.0)
        x[n + 1] = 0.0; y[n + 1] =   a; z[n + 1] =   a; w[n + 1] = v;
        x[n + 2] = 0.0; y[n + 2] =   a; z[n + 2] =  -a; w[n + 2] = v;
        x[n + 3] = 0.0; y[n + 3] =  -a; z[n + 3] =   a; w[n + 3] = v;
        x[n + 4] = 0.0; y[n + 4] =  -a; z[n + 4] =  -a; w[n + 4] = v;
        x[n + 5] =   a; y[n + 5] = 0.0; z[n + 5] =   a; w[n + 5] = v;
        x[n + 6] =   a; y[n + 6] = 0.0; z[n + 6] =  -a; w[n + 6] = v;
        x[n + 7] =  -a; y[n + 7] = 0.0; z[n + 7] =   a; w[n + 7] = v;
        x[n + 8] =  -a; y[n + 8] = 0.0; z[n + 8] =  -a; w[n + 8] = v;
        x[n + 9] =   a; y[n + 9] =   a; z[n + 9] = 0.0; w[n + 9] = v;
        x[n + 10] =   a; y[n + 10] =  -a; z[n + 10] = 0.0; w[n + 10] = v;
        x[n + 11] =  -a; y[n + 11] =   a; z[n + 11] = 0.0; w[n + 11] = v;
        x[n + 12] =  -a; y[n + 12] =  -a; z[n + 12] = 0.0; w[n + 12] = v;
        n = 12
    elseif code == 3
        a = 1.0 / sqrt(3.0)
        x[n + 1] =   a; y[n + 1] =   a; z[n + 1] =   a; w[n + 1] = v;
        x[n + 2] =   a; y[n + 2] =   a; z[n + 2] =  -a; w[n + 2] = v;
        x[n + 3] =   a; y[n + 3] =  -a; z[n + 3] =   a; w[n + 3] = v;
        x[n + 4] =   a; y[n + 4] =  -a; z[n + 4] =  -a; w[n + 4] = v;
        x[n + 5] =  -a; y[n + 5] =   a; z[n + 5] =   a; w[n + 5] = v;
        x[n + 6] =  -a; y[n + 6] =   a; z[n + 6] =  -a; w[n + 6] = v;
        x[n + 7] =  -a; y[n + 7] =  -a; z[n + 7] =   a; w[n + 7] = v;
        x[n + 8] =  -a; y[n + 8] =  -a; z[n + 8] =  -a; w[n + 8] = v;
        n = 8
    elseif code == 4
        b = sqrt(1.0 - 2.0 * a^2)
        x[n + 1] =   a; y[n + 1] =   a; z[n + 1] =   b; w[n + 1] = v;
        x[n + 2] =   a; y[n + 2] =   a; z[n + 2] =  -b; w[n + 2] = v;
        x[n + 3] =   a; y[n + 3] =  -a; z[n + 3] =   b; w[n + 3] = v;
        x[n + 4] =   a; y[n + 4] =  -a; z[n + 4] =  -b; w[n + 4] = v;
        x[n + 5] =  -a; y[n + 5] =   a; z[n + 5] =   b; w[n + 5] = v;
        x[n + 6] =  -a; y[n + 6] =   a; z[n + 6] =  -b; w[n + 6] = v;
        x[n + 7] =  -a; y[n + 7] =  -a; z[n + 7] =   b; w[n + 7] = v;
        x[n + 8] =  -a; y[n + 8] =  -a; z[n + 8] =  -b; w[n + 8] = v;
        x[n + 9] =   a; y[n + 9] =   b; z[n + 9] =   a; w[n + 9] = v;
        x[n + 10] =   a; y[n + 10] =  -b; z[n + 10] =   a; w[n + 10] = v;
        x[n + 11] =   a; y[n + 11] =   b; z[n + 11] =  -a; w[n + 11] = v;
        x[n + 12] =   a; y[n + 12] =  -b; z[n + 12] =  -a; w[n + 12] = v;
        x[n + 13] =  -a; y[n + 13] =   b; z[n + 13] =   a; w[n + 13] = v;
        x[n + 14] =  -a; y[n + 14] =  -b; z[n + 14] =   a; w[n + 14] = v;
        x[n + 15] =  -a; y[n + 15] =   b; z[n + 15] =  -a; w[n + 15] = v;
        x[n + 16] =  -a; y[n + 16] =  -b; z[n + 16] =  -a; w[n + 16] = v;
        x[n + 17] =   b; y[n + 17] =   a; z[n + 17] =   a; w[n + 17] = v;
        x[n + 18] =  -b; y[n + 18] =   a; z[n + 18] =   a; w[n + 18] = v;
        x[n + 19] =   b; y[n + 19] =   a; z[n + 19] =  -a; w[n + 19] = v;
        x[n + 20] =  -b; y[n + 20] =   a; z[n + 20] =  -a; w[n + 20] = v;
        x[n + 21] =   b; y[n + 21] =  -a; z[n + 21] =   a; w[n + 21] = v;
        x[n + 22] =  -b; y[n + 22] =  -a; z[n + 22] =   a; w[n + 22] = v;
        x[n + 23] =   b; y[n + 23] =  -a; z[n + 23] =  -a; w[n + 23] = v;
        x[n + 24] =  -b; y[n + 24] =  -a; z[n + 24] =  -a; w[n + 24] = v;
        n = 24
    elseif code == 5
        b = sqrt(1.0 - a^2)
        x[n + 1] =   a; y[n + 1] =   b; z[n + 1] = 0.0; w[n + 1] = v;
        x[n + 2] =   a; y[n + 2] =  -b; z[n + 2] = 0.0; w[n + 2] = v;
        x[n + 3] =  -a; y[n + 3] =   b; z[n + 3] = 0.0; w[n + 3] = v;
        x[n + 4] =  -a; y[n + 4] =  -b; z[n + 4] = 0.0; w[n + 4] = v;
        x[n + 5] =   b; y[n + 5] =   a; z[n + 5] = 0.0; w[n + 5] = v;
        x[n + 6] =   b; y[n + 6] =  -a; z[n + 6] = 0.0; w[n + 6] = v;
        x[n + 7] =  -b; y[n + 7] =   a; z[n + 7] = 0.0; w[n + 7] = v;
        x[n + 8] =  -b; y[n + 8] =  -a; z[n + 8] = 0.0; w[n + 8] = v;
        x[n + 9] =   a; y[n + 9] = 0.0; z[n + 9] =   b; w[n + 9] = v;
        x[n + 10] =   a; y[n + 10] = 0.0; z[n + 10] =  -b; w[n + 10] = v;
        x[n + 11] =  -a; y[n + 11] = 0.0; z[n + 11] =   b; w[n + 11] = v;
        x[n + 12] =  -a; y[n + 12] = 0.0; z[n + 12] =  -b; w[n + 12] = v;
        x[n + 13] =   b; y[n + 13] = 0.0; z[n + 13] =   a; w[n + 13] = v;
        x[n + 14] =   b; y[n + 14] = 0.0; z[n + 14] =  -a; w[n + 14] = v;
        x[n + 15] =  -b; y[n + 15] = 0.0; z[n + 15] =   a; w[n + 15] = v;
        x[n + 16] =  -b; y[n + 16] = 0.0; z[n + 16] =  -a; w[n + 16] = v;
        x[n + 17] = 0.0; y[n + 17] =   a; z[n + 17] =   b; w[n + 17] = v;
        x[n + 18] = 0.0; y[n + 18] =   a; z[n + 18] =  -b; w[n + 18] = v;
        x[n + 19] = 0.0; y[n + 19] =  -a; z[n + 19] =   b; w[n + 19] = v;
        x[n + 20] = 0.0; y[n + 20] =  -a; z[n + 20] =  -b; w[n + 20] = v;
        x[n + 21] = 0.0; y[n + 21] =   b; z[n + 21] =   a; w[n + 21] = v;
        x[n + 22] = 0.0; y[n + 22] =   b; z[n + 22] =  -a; w[n + 22] = v;
        x[n + 23] = 0.0; y[n + 23] =  -b; z[n + 23] =   a; w[n + 23] = v;
        x[n + 24] = 0.0; y[n + 24] =  -b; z[n + 24] =  -a; w[n + 24] = v;
        n = 24
    elseif code == 6
        c = sqrt(1.0 - a^2 - b^2)
        x[n + 1] =   a; y[n + 1] =   b; z[n + 1] =   c; w[n + 1] = v;
        x[n + 2] =   a; y[n + 2] =   b; z[n + 2] =  -c; w[n + 2] = v;
        x[n + 3] =   a; y[n + 3] =  -b; z[n + 3] =   c; w[n + 3] = v;
        x[n + 4] =   a; y[n + 4] =  -b; z[n + 4] =  -c; w[n + 4] = v;
        x[n + 5] =  -a; y[n + 5] =   b; z[n + 5] =   c; w[n + 5] = v;
        x[n + 6] =  -a; y[n + 6] =   b; z[n + 6] =  -c; w[n + 6] = v;
        x[n + 7] =  -a; y[n + 7] =  -b; z[n + 7] =   c; w[n + 7] = v;
        x[n + 8] =  -a; y[n + 8] =  -b; z[n + 8] =  -c; w[n + 8] = v;
        x[n + 9] =   a; y[n + 9] =   c; z[n + 9] =   b; w[n + 9] = v;
        x[n + 10] =   a; y[n + 10] =   c; z[n + 10] =  -b; w[n + 10] = v;
        x[n + 11] =   a; y[n + 11] =  -c; z[n + 11] =   b; w[n + 11] = v;
        x[n + 12] =   a; y[n + 12] =  -c; z[n + 12] =  -b; w[n + 12] = v;
        x[n + 13] =  -a; y[n + 13] =   c; z[n + 13] =   b; w[n + 13] = v;
        x[n + 14] =  -a; y[n + 14] =   c; z[n + 14] =  -b; w[n + 14] = v;
        x[n + 15] =  -a; y[n + 15] =  -c; z[n + 15] =   b; w[n + 15] = v;
        x[n + 16] =  -a; y[n + 16] =  -c; z[n + 16] =  -b; w[n + 16] = v;
        x[n + 17] =   b; y[n + 17] =   a; z[n + 17] =   c; w[n + 17] = v;
        x[n + 18] =   b; y[n + 18] =   a; z[n + 18] =  -c; w[n + 18] = v;
        x[n + 19] =   b; y[n + 19] =  -a; z[n + 19] =   c; w[n + 19] = v;
        x[n + 20] =   b; y[n + 20] =  -a; z[n + 20] =  -c; w[n + 20] = v;
        x[n + 21] =  -b; y[n + 21] =   a; z[n + 21] =   c; w[n + 21] = v;
        x[n + 22] =  -b; y[n + 22] =   a; z[n + 22] =  -c; w[n + 22] = v;
        x[n + 23] =  -b; y[n + 23] =  -a; z[n + 23] =   c; w[n + 23] = v;
        x[n + 24] =  -b; y[n + 24] =  -a; z[n + 24] =  -c; w[n + 24] = v;
        x[n + 25] =   b; y[n + 25] =   c; z[n + 25] =   a; w[n + 25] = v;
        x[n + 26] =   b; y[n + 26] =   c; z[n + 26] =  -a; w[n + 26] = v;
        x[n + 27] =   b; y[n + 27] =  -c; z[n + 27] =   a; w[n + 27] = v;
        x[n + 28] =   b; y[n + 28] =  -c; z[n + 28] =  -a; w[n + 28] = v;
        x[n + 29] =  -b; y[n + 29] =   c; z[n + 29] =   a; w[n + 29] = v;
        x[n + 30] =  -b; y[n + 30] =   c; z[n + 30] =  -a; w[n + 30] = v;
        x[n + 31] =  -b; y[n + 31] =  -c; z[n + 31] =   a; w[n + 31] = v;
        x[n + 32] =  -b; y[n + 32] =  -c; z[n + 32] =  -a; w[n + 32] = v;
        x[n + 33] =   c; y[n + 33] =   a; z[n + 33] =   b; w[n + 33] = v;
        x[n + 34] =   c; y[n + 34] =   a; z[n + 34] =  -b; w[n + 34] = v;
        x[n + 35] =   c; y[n + 35] =  -a; z[n + 35] =   b; w[n + 35] = v;
        x[n + 36] =   c; y[n + 36] =  -a; z[n + 36] =  -b; w[n + 36] = v;
        x[n + 37] =  -c; y[n + 37] =   a; z[n + 37] =   b; w[n + 37] = v;
        x[n + 38] =  -c; y[n + 38] =   a; z[n + 38] =  -b; w[n + 38] = v;
        x[n + 39] =  -c; y[n + 39] =  -a; z[n + 39] =   b; w[n + 39] = v;
        x[n + 40] =  -c; y[n + 40] =  -a; z[n + 40] =  -b; w[n + 40] = v;
        x[n + 41] =   c; y[n + 41] =   b; z[n + 41] =   a; w[n + 41] = v;
        x[n + 42] =   c; y[n + 42] =   b; z[n + 42] =  -a; w[n + 42] = v;
        x[n + 43] =   c; y[n + 43] =  -b; z[n + 43] =   a; w[n + 43] = v;
        x[n + 44] =   c; y[n + 44] =  -b; z[n + 44] =  -a; w[n + 44] = v;
        x[n + 45] =  -c; y[n + 45] =   b; z[n + 45] =   a; w[n + 45] = v;
        x[n + 46] =  -c; y[n + 46] =   b; z[n + 46] =  -a; w[n + 46] = v;
        x[n + 47] =  -c; y[n + 47] =  -b; z[n + 47] =   a; w[n + 47] = v;
        x[n + 48] =  -c; y[n + 48] =  -b; z[n + 48] =  -a; w[n + 48] = v;
        n = 48
    else
        throw(ArgumentError("code must be 1, 2, 3, 4, 5 or 6"))
    end

    return n
end


"""
    ld0006!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 6 point Lebedev angular grid.
"""
function ld0006!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.1666666666666667
    n += gen_oh!(1, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0014!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 14 point Lebedev angular grid.
"""
function ld0014!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.6666666666666667e-1
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.7500000000000000e-1
    n += gen_oh!(3, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0026!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 26 point Lebedev angular grid.
"""
function ld0026!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.4761904761904762e-1
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.3809523809523810e-1
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.3214285714285714e-1
    n += gen_oh!(3, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0038!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 38 point Lebedev angular grid.
"""
function ld0038!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.9523809523809524e-2
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.3214285714285714e-1
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.4597008433809831
    v  = 0.2857142857142857e-1
    n += gen_oh!(5, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0050!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 50 point Lebedev angular grid.
"""
function ld0050!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.1269841269841270e-1
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.2257495590828924e-1
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.2109375000000000e-1
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.3015113445777636
    v  = 0.2017333553791887e-1
    n += gen_oh!(4, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0074!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 74 point Lebedev angular grid.
"""
function ld0074!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.5130671797338464e-3
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.1660406956574204e-1
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = -0.2958603896103896e-1
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.4803844614152614
    v  = 0.2657620708215946e-1
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3207726489807764
    v  = 0.1652217099371571e-1
    n += gen_oh!(5, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0086!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 86 point Lebedev angular grid.
"""
function ld0086!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.1154401154401154e-1
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.1194390908585628e-1
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.3696028464541502
    v  = 0.1111055571060340e-1
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6943540066026664
    v  = 0.1187650129453714e-1
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3742430390903412
    v  = 0.1181230374690448e-1
    n += gen_oh!(5, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0110!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 110 point Lebedev angular grid.
"""
function ld0110!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.3828270494937162e-2
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.9793737512487512e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.1851156353447362
    v  = 0.8211737283191111e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6904210483822922
    v  = 0.9942814891178103e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3956894730559419
    v  = 0.9595471336070963e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4783690288121502
    v  = 0.9694996361663028e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0146!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 146 point Lebedev angular grid.
"""
function ld0146!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.5996313688621381e-3
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.7372999718620756e-2
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.7210515360144488e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.6764410400114264
    v  = 0.7116355493117555e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4174961227965453
    v  = 0.6753829486314477e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1574676672039082
    v  = 0.7574394159054034e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1403553811713183
    b  = 0.4493328323269557
    v  = 0.6991087353303262e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0170!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 170 point Lebedev angular grid.
"""
function ld0170!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.5544842902037365e-2
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.6071332770670752e-2
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.6383674773515093e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.2551252621114134
    v  = 0.5183387587747790e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6743601460362766
    v  = 0.6317929009813725e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4318910696719410
    v  = 0.6201670006589077e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2613931360335988
    v  = 0.5477143385137348e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.4990453161796037
    b  = 0.1446630744325115
    v  = 0.5968383987681156e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0194!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 194 point Lebedev angular grid.
"""
function ld0194!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.1782340447244611e-2
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.5716905949977102e-2
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.5573383178848738e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.6712973442695226
    v  = 0.5608704082587997e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2892465627575439
    v  = 0.5158237711805383e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4446933178717437
    v  = 0.5518771467273614e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1299335447650067
    v  = 0.4106777028169394e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3457702197611283
    v  = 0.5051846064614808e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1590417105383530
    b  = 0.8360360154824589
    v  = 0.5530248916233094e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0230!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 230 point Lebedev angular grid.
"""
function ld0230!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = -0.5522639919727325e-1
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.4450274607445226e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.4492044687397611
    v  = 0.4496841067921404e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2520419490210201
    v  = 0.5049153450478750e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6981906658447242
    v  = 0.3976408018051883e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6587405243460960
    v  = 0.4401400650381014e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4038544050097660e-1
    v  = 0.1724544350544401e-1
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5823842309715585
    v  = 0.4231083095357343e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3545877390518688
    v  = 0.5198069864064399e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2272181808998187
    b  = 0.4864661535886647
    v  = 0.4695720972568883e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0266!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 266 point Lebedev angular grid.
"""
function ld0266!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = -0.1313769127326952e-2
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = -0.2522728704859336e-2
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.4186853881700583e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.7039373391585475
    v  = 0.5315167977810885e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1012526248572414
    v  = 0.4047142377086219e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4647448726420539
    v  = 0.4112482394406990e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3277420654971629
    v  = 0.3595584899758782e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6620338663699974
    v  = 0.4256131351428158e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.8506508083520399
    v  = 0.4229582700647240e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3233484542692899
    b  = 0.1153112011009701
    v  = 0.4080914225780505e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2314790158712601
    b  = 0.5244939240922365
    v  = 0.4071467593830964e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0302!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 302 point Lebedev angular grid.
"""
function ld0302!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.8545911725128148e-3
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.3599119285025571e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.3515640345570105
    v  = 0.3449788424305883e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6566329410219612
    v  = 0.3604822601419882e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4729054132581005
    v  = 0.3576729661743367e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.9618308522614784e-1
    v  = 0.2352101413689164e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2219645236294178
    v  = 0.3108953122413675e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7011766416089545
    v  = 0.3650045807677255e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2644152887060663
    v  = 0.2982344963171804e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5718955891878961
    v  = 0.3600820932216460e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2510034751770465
    b  = 0.8000727494073952
    v  = 0.3571540554273387e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1233548532583327
    b  = 0.4127724083168531
    v  = 0.3392312205006170e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0350!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 350 point Lebedev angular grid.
"""
function ld0350!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.3006796749453936e-2
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.3050627745650771e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.7068965463912316
    v  = 0.1621104600288991e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4794682625712025
    v  = 0.3005701484901752e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1927533154878019
    v  = 0.2990992529653774e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6930357961327123
    v  = 0.2982170644107595e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3608302115520091
    v  = 0.2721564237310992e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6498486161496169
    v  = 0.3033513795811141e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1932945013230339
    v  = 0.3007949555218533e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3800494919899303
    v  = 0.2881964603055307e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2899558825499574
    b  = 0.7934537856582316
    v  = 0.2958357626535696e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.9684121455103957e-1
    b  = 0.8280801506686862
    v  = 0.3036020026407088e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1833434647041659
    b  = 0.9074658265305127
    v  = 0.2832187403926303e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0434!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 434 point Lebedev angular grid.
"""
function ld0434!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.5265897968224436e-3
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.2548219972002607e-2
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.2512317418927307e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.6909346307509111
    v  = 0.2530403801186355e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1774836054609158
    v  = 0.2014279020918528e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4914342637784746
    v  = 0.2501725168402936e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6456664707424256
    v  = 0.2513267174597564e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2861289010307638
    v  = 0.2302694782227416e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7568084367178018e-1
    v  = 0.1462495621594614e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3927259763368002
    v  = 0.2445373437312980e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.8818132877794288
    v  = 0.2417442375638981e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.9776428111182649
    v  = 0.1910951282179532e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2054823696403044
    b  = 0.8689460322872412
    v  = 0.2416930044324775e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5905157048925271
    b  = 0.7999278543857286
    v  = 0.2512236854563495e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5550152361076807
    b  = 0.7717462626915901
    v  = 0.2496644054553086e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.9371809858553722
    b  = 0.3344363145343455
    v  = 0.2236607760437849e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0590!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 590 point Lebedev angular grid.
"""
function ld0590!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.3095121295306187e-3
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.1852379698597489e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.7040954938227469
    v  = 0.1871790639277744e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6807744066455243
    v  = 0.1858812585438317e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6372546939258752
    v  = 0.1852028828296213e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5044419707800358
    v  = 0.1846715956151242e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4215761784010967
    v  = 0.1818471778162769e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3317920736472123
    v  = 0.1749564657281154e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2384736701421887
    v  = 0.1617210647254411e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1459036449157763
    v  = 0.1384737234851692e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6095034115507196e-1
    v  = 0.9764331165051050e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6116843442009876
    v  = 0.1857161196774078e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3964755348199858
    v  = 0.1705153996395864e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1724782009907724
    v  = 0.1300321685886048e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5610263808622060
    b  = 0.3518280927733519
    v  = 0.1842866472905286e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4742392842551980
    b  = 0.2634716655937950
    v  = 0.1802658934377451e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5984126497885380
    b  = 0.1816640840360209
    v  = 0.1849830560443660e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3791035407695563
    b  = 0.1720795225656878
    v  = 0.1713904507106709e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2778673190586244
    b  = 0.8213021581932511e-1
    v  = 0.1555213603396808e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5033564271075117
    b  = 0.8999205842074875e-1
    v  = 0.1802239128008525e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0770!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 770 point Lebedev angular grid.
"""
function ld0770!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.2192942088181184e-3
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.1436433617319080e-2
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.1421940344335877e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.5087204410502360e-1
    v  = 0.6798123511050502e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1228198790178831
    v  = 0.9913184235294912e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2026890814408786
    v  = 0.1180207833238949e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2847745156464294
    v  = 0.1296599602080921e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3656719078978026
    v  = 0.1365871427428316e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4428264886713469
    v  = 0.1402988604775325e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5140619627249735
    v  = 0.1418645563595609e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6306401219166803
    v  = 0.1421376741851662e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6716883332022612
    v  = 0.1423996475490962e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6979792685336881
    v  = 0.1431554042178567e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1446865674195309
    v  = 0.9254401499865368e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3390263475411216
    v  = 0.1250239995053509e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5335804651263506
    v  = 0.1394365843329230e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.6944024393349413e-1
    b  = 0.2355187894242326
    v  = 0.1127089094671749e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2269004109529460
    b  = 0.4102182474045730
    v  = 0.1345753760910670e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.8025574607775339e-1
    b  = 0.6214302417481605
    v  = 0.1424957283316783e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1467999527896572
    b  = 0.3245284345717394
    v  = 0.1261523341237750e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1571507769824727
    b  = 0.5224482189696630
    v  = 0.1392547106052696e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2365702993157246
    b  = 0.6017546634089558
    v  = 0.1418761677877656e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.7714815866765732e-1
    b  = 0.4346575516141163
    v  = 0.1338366684479554e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3062936666210730
    b  = 0.4908826589037616
    v  = 0.1393700862676131e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3822477379524787
    b  = 0.5648768149099500
    v  = 0.1415914757466932e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld0974!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 974 point Lebedev angular grid.
"""
function ld0974!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.1438294190527431e-3
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.1125772288287004e-2
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.4292963545341347e-1
    v  = 0.4948029341949241e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1051426854086404
    v  = 0.7357990109125470e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1750024867623087
    v  = 0.8889132771304384e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2477653379650257
    v  = 0.9888347838921435e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3206567123955957
    v  = 0.1053299681709471e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3916520749849983
    v  = 0.1092778807014578e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4590825874187624
    v  = 0.1114389394063227e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5214563888415861
    v  = 0.1123724788051555e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6253170244654199
    v  = 0.1125239325243814e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6637926744523170
    v  = 0.1126153271815905e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6910410398498301
    v  = 0.1130286931123841e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7052907007457760
    v  = 0.1134986534363955e-2
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1236686762657990
    v  = 0.6823367927109931e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2940777114468387
    v  = 0.9454158160447096e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.4697753849207649
    v  = 0.1074429975385679e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.6334563241139567
    v  = 0.1129300086569132e-2
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5974048614181342e-1
    b  = 0.2029128752777523
    v  = 0.8436884500901954e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1375760408473636
    b  = 0.4602621942484054
    v  = 0.1075255720448885e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3391016526336286
    b  = 0.5030673999662036
    v  = 0.1108577236864462e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1271675191439820
    b  = 0.2817606422442134
    v  = 0.9566475323783357e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2693120740413512
    b  = 0.4331561291720157
    v  = 0.1080663250717391e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1419786452601918
    b  = 0.6256167358580814
    v  = 0.1126797131196295e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6709284600738255e-1
    b  = 0.3798395216859157
    v  = 0.1022568715358061e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.7057738183256172e-1
    b  = 0.5517505421423520
    v  = 0.1108960267713108e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2783888477882155
    b  = 0.6029619156159187
    v  = 0.1122790653435766e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1979578938917407
    b  = 0.3589606329589096
    v  = 0.1032401847117460e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2087307061103274
    b  = 0.5348666438135476
    v  = 0.1107249382283854e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4055122137872836
    b  = 0.5674997546074373
    v  = 0.1121780048519972e-2
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld1202!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 1202 point Lebedev angular grid.
"""
function ld1202!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.1105189233267572e-3
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.9205232738090741e-3
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.9133159786443561e-3
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.3712636449657089e-1
    v  = 0.3690421898017899e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.9140060412262223e-1
    v  = 0.5603990928680660e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1531077852469906
    v  = 0.6865297629282609e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2180928891660612
    v  = 0.7720338551145630e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2839874532200175
    v  = 0.8301545958894795e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3491177600963764
    v  = 0.8686692550179628e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4121431461444309
    v  = 0.8927076285846890e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4718993627149127
    v  = 0.9060820238568219e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5273145452842337
    v  = 0.9119777254940867e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6209475332444019
    v  = 0.9128720138604181e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6569722711857291
    v  = 0.9130714935691735e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6841788309070143
    v  = 0.9152873784554116e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7012604330123631
    v  = 0.9187436274321654e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1072382215478166
    v  = 0.5176977312965694e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2582068959496968
    v  = 0.7331143682101417e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.4172752955306717
    v  = 0.8463232836379928e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5700366911792503
    v  = 0.9031122694253992e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.9827986018263947
    b  = 0.1771774022615325
    v  = 0.6485778453163257e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.9624249230326228
    b  = 0.2475716463426288
    v  = 0.7435030910982369e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.9402007994128811
    b  = 0.3354616289066489
    v  = 0.7998527891839054e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.9320822040143202
    b  = 0.3173615246611977
    v  = 0.8101731497468018e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.9043674199393299
    b  = 0.4090268427085357
    v  = 0.8483389574594331e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.8912407560074747
    b  = 0.3854291150669224
    v  = 0.8556299257311812e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.8676435628462708
    b  = 0.4932221184851285
    v  = 0.8803208679738260e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.8581979986041619
    b  = 0.4785320675922435
    v  = 0.8811048182425720e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.8396753624049856
    b  = 0.4507422593157064
    v  = 0.8850282341265444e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.8165288564022188
    b  = 0.5632123020762100
    v  = 0.9021342299040653e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.8015469370783529
    b  = 0.5434303569693900
    v  = 0.9010091677105086e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.7773563069070351
    b  = 0.5123518486419871
    v  = 0.9022692938426915e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.7661621213900394
    b  = 0.6394279634749102
    v  = 0.9158016174693465e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.7553584143533510
    b  = 0.6269805509024392
    v  = 0.9131578003189435e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.7344305757559503
    b  = 0.6031161693096310
    v  = 0.9107813579482705e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.7043837184021765
    b  = 0.5693702498468441
    v  = 0.9105760258970126e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld1454!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 1454 point Lebedev angular grid.
"""
function ld1454!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.7777160743261247e-4
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.7557646413004701e-3
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.3229290663413854e-1
    v  = 0.2841633806090617e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.8036733271462222e-1
    v  = 0.4374419127053555e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1354289960531653
    v  = 0.5417174740872172e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1938963861114426
    v  = 0.6148000891358593e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2537343715011275
    v  = 0.6664394485800705e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3135251434752570
    v  = 0.7025039356923220e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3721558339375338
    v  = 0.7268511789249627e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4286809575195696
    v  = 0.7422637534208629e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4822510128282994
    v  = 0.7509545035841214e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5320679333566263
    v  = 0.7548535057718401e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6172998195394274
    v  = 0.7554088969774001e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6510679849127481
    v  = 0.7553147174442808e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6777315251687360
    v  = 0.7564767653292297e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6963109410648741
    v  = 0.7587991808518730e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7058935009831749
    v  = 0.7608261832033027e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.9955546194091857
    v  = 0.4021680447874916e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.9734115901794209
    v  = 0.5804871793945964e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.9275693732388626
    v  = 0.6792151955945159e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.8568022422795103
    v  = 0.7336741211286294e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.7623495553719372
    v  = 0.7581866300989608e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5707522908892223
    b  = 0.4387028039889501
    v  = 0.7538257859800743e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5196463388403083
    b  = 0.3858908414762617
    v  = 0.7483517247053123e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4646337531215351
    b  = 0.3301937372343854
    v  = 0.7371763661112059e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4063901697557691
    b  = 0.2725423573563777
    v  = 0.7183448895756934e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3456329466643087
    b  = 0.2139510237495250
    v  = 0.6895815529822191e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2831395121050332
    b  = 0.1555922309786647
    v  = 0.6480105801792886e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2197682022925330
    b  = 0.9892878979686097e-1
    v  = 0.5897558896594636e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1564696098650355
    b  = 0.4598642910675510e-1
    v  = 0.5095708849247346e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6027356673721295
    b  = 0.3376625140173426
    v  = 0.7536906428909755e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5496032320255096
    b  = 0.2822301309727988
    v  = 0.7472505965575118e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4921707755234567
    b  = 0.2248632342592540
    v  = 0.7343017132279698e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4309422998598483
    b  = 0.1666224723456479
    v  = 0.7130871582177445e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3664108182313672
    b  = 0.1086964901822169
    v  = 0.6817022032112776e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2990189057758436
    b  = 0.5251989784120085e-1
    v  = 0.6380941145604121e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6268724013144998
    b  = 0.2297523657550023
    v  = 0.7550381377920310e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5707324144834607
    b  = 0.1723080607093800
    v  = 0.7478646640144802e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5096360901960365
    b  = 0.1140238465390513
    v  = 0.7335918720601220e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4438729938312456
    b  = 0.5611522095882537e-1
    v  = 0.7110120527658118e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6419978471082389
    b  = 0.1164174423140873
    v  = 0.7571363978689501e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5817218061802611
    b  = 0.5797589531445219e-1
    v  = 0.7489908329079234e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld1730!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 1730 point Lebedev angular grid.
"""
function ld1730!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.6309049437420976e-4
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.6398287705571748e-3
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.6357185073530720e-3
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.2860923126194662e-1
    v  = 0.2221207162188168e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7142556767711522e-1
    v  = 0.3475784022286848e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1209199540995559
    v  = 0.4350742443589804e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1738673106594379
    v  = 0.4978569136522127e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2284645438467734
    v  = 0.5435036221998053e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2834807671701512
    v  = 0.5765913388219542e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3379680145467339
    v  = 0.6001200359226003e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3911355454819537
    v  = 0.6162178172717512e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4422860353001403
    v  = 0.6265218152438485e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4907781568726057
    v  = 0.6323987160974212e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5360006153211468
    v  = 0.6350767851540569e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6142105973596603
    v  = 0.6354362775297107e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6459300387977504
    v  = 0.6352302462706235e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6718056125089225
    v  = 0.6358117881417972e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6910888533186254
    v  = 0.6373101590310117e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7030467416823252
    v  = 0.6390428961368665e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.8354951166354646e-1
    v  = 0.3186913449946576e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2050143009099486
    v  = 0.4678028558591711e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3370208290706637
    v  = 0.5538829697598626e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.4689051484233963
    v  = 0.6044475907190476e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5939400424557334
    v  = 0.6313575103509012e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1394983311832261
    b  = 0.4097581162050343e-1
    v  = 0.4078626431855630e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1967999180485014
    b  = 0.8851987391293348e-1
    v  = 0.4759933057812725e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2546183732548967
    b  = 0.1397680182969819
    v  = 0.5268151186413440e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3121281074713875
    b  = 0.1929452542226526
    v  = 0.5643048560507316e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3685981078502492
    b  = 0.2467898337061562
    v  = 0.5914501076613073e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4233760321547856
    b  = 0.3003104124785409
    v  = 0.6104561257874195e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4758671236059246
    b  = 0.3526684328175033
    v  = 0.6230252860707806e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5255178579796463
    b  = 0.4031134861145713
    v  = 0.6305618761760796e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5718025633734589
    b  = 0.4509426448342351
    v  = 0.6343092767597889e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2686927772723415
    b  = 0.4711322502423248e-1
    v  = 0.5176268945737826e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3306006819904809
    b  = 0.9784487303942695e-1
    v  = 0.5564840313313692e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3904906850594983
    b  = 0.1505395810025273
    v  = 0.5856426671038980e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4479957951904390
    b  = 0.2039728156296050
    v  = 0.6066386925777091e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5027076848919780
    b  = 0.2571529941121107
    v  = 0.6208824962234458e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5542087392260217
    b  = 0.3092191375815670
    v  = 0.6296314297822907e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6020850887375187
    b  = 0.3593807506130276
    v  = 0.6340423756791859e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4019851409179594
    b  = 0.5063389934378671e-1
    v  = 0.5829627677107342e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4635614567449800
    b  = 0.1032422269160612
    v  = 0.6048693376081110e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5215860931591575
    b  = 0.1566322094006254
    v  = 0.6202362317732461e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5758202499099271
    b  = 0.2098082827491099
    v  = 0.6299005328403779e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6259893683876795
    b  = 0.2618824114553391
    v  = 0.6347722390609353e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5313795124811891
    b  = 0.5263245019338556e-1
    v  = 0.6203778981238834e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5893317955931995
    b  = 0.1061059730982005
    v  = 0.6308414671239979e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6426246321215801
    b  = 0.1594171564034221
    v  = 0.6362706466959498e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6511904367376113
    b  = 0.5354789536565540e-1
    v  = 0.6375414170333233e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld2030!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 2030 point Lebedev angular grid.
"""
function ld2030!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.4656031899197431e-4
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.5421549195295507e-3
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.2540835336814348e-1
    v  = 0.1778522133346553e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6399322800504915e-1
    v  = 0.2811325405682796e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1088269469804125
    v  = 0.3548896312631459e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1570670798818287
    v  = 0.4090310897173364e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2071163932282514
    v  = 0.4493286134169965e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2578914044450844
    v  = 0.4793728447962723e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3085687558169623
    v  = 0.5015415319164265e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3584719706267024
    v  = 0.5175127372677937e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4070135594428709
    v  = 0.5285522262081019e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4536618626222638
    v  = 0.5356832703713962e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4979195686463577
    v  = 0.5397914736175170e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5393075111126999
    v  = 0.5416899441599930e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6115617676843916
    v  = 0.5419308476889938e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6414308435160159
    v  = 0.5416936902030596e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6664099412721607
    v  = 0.5419544338703164e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6859161771214913
    v  = 0.5428983656630975e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6993625593503890
    v  = 0.5442286500098193e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7062393387719380
    v  = 0.5452250345057301e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7479028168349763e-1
    v  = 0.2568002497728530e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1848951153969366
    v  = 0.3827211700292145e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3059529066581305
    v  = 0.4579491561917824e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.4285556101021362
    v  = 0.5042003969083574e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5468758653496526
    v  = 0.5312708889976025e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.6565821978343439
    v  = 0.5438401790747117e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1253901572367117
    b  = 0.3681917226439641e-1
    v  = 0.3316041873197344e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1775721510383941
    b  = 0.7982487607213301e-1
    v  = 0.3899113567153771e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2305693358216114
    b  = 0.1264640966592335
    v  = 0.4343343327201309e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2836502845992063
    b  = 0.1751585683418957
    v  = 0.4679415262318919e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3361794746232590
    b  = 0.2247995907632670
    v  = 0.4930847981631031e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3875979172264824
    b  = 0.2745299257422246
    v  = 0.5115031867540091e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4374019316999074
    b  = 0.3236373482441118
    v  = 0.5245217148457367e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4851275843340022
    b  = 0.3714967859436741
    v  = 0.5332041499895321e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5303391803806868
    b  = 0.4175353646321745
    v  = 0.5384583126021542e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5726197380596287
    b  = 0.4612084406355461
    v  = 0.5411067210798852e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2431520732564863
    b  = 0.4258040133043952e-1
    v  = 0.4259797391468714e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3002096800895869
    b  = 0.8869424306722721e-1
    v  = 0.4604931368460021e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3558554457457432
    b  = 0.1368811706510655
    v  = 0.4871814878255202e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4097782537048887
    b  = 0.1860739985015033
    v  = 0.5072242910074885e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4616337666067458
    b  = 0.2354235077395853
    v  = 0.5217069845235350e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5110707008417874
    b  = 0.2842074921347011
    v  = 0.5315785966280310e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5577415286163795
    b  = 0.3317784414984102
    v  = 0.5376833708758905e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6013060431366950
    b  = 0.3775299002040700
    v  = 0.5408032092069521e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3661596767261781
    b  = 0.4599367887164592e-1
    v  = 0.4842744917904866e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4237633153506581
    b  = 0.9404893773654421e-1
    v  = 0.5048926076188130e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4786328454658452
    b  = 0.1431377109091971
    v  = 0.5202607980478373e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5305702076789774
    b  = 0.1924186388843570
    v  = 0.5309932388325743e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5793436224231788
    b  = 0.2411590944775190
    v  = 0.5377419770895208e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6247069017094747
    b  = 0.2886871491583605
    v  = 0.5411696331677717e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4874315552535204
    b  = 0.4804978774953206e-1
    v  = 0.5197996293282420e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5427337322059053
    b  = 0.9716857199366665e-1
    v  = 0.5311120836622945e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5943493747246700
    b  = 0.1465205839795055
    v  = 0.5384309319956951e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6421314033564943
    b  = 0.1953579449803574
    v  = 0.5421859504051886e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6020628374713980
    b  = 0.4916375015738108e-1
    v  = 0.5390948355046314e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6529222529856881
    b  = 0.9861621540127005e-1
    v  = 0.5433312705027845e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld2354!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 2354 point Lebedev angular grid.
"""
function ld2354!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.3922616270665292e-4
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.4703831750854424e-3
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.4678202801282136e-3
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.2290024646530589e-1
    v  = 0.1437832228979900e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5779086652271284e-1
    v  = 0.2303572493577644e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.9863103576375984e-1
    v  = 0.2933110752447454e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1428155792982185
    v  = 0.3402905998359838e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1888978116601463
    v  = 0.3759138466870372e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2359091682970210
    v  = 0.4030638447899798e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2831228833706171
    v  = 0.4236591432242211e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3299495857966693
    v  = 0.4390522656946746e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3758840802660796
    v  = 0.4502523466626247e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4204751831009480
    v  = 0.4580577727783541e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4633068518751051
    v  = 0.4631391616615899e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5039849474507313
    v  = 0.4660928953698676e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5421265793440747
    v  = 0.4674751807936953e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6092660230557310
    v  = 0.4676414903932920e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6374654204984869
    v  = 0.4674086492347870e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6615136472609892
    v  = 0.4674928539483207e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6809487285958127
    v  = 0.4680748979686447e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6952980021665196
    v  = 0.4690449806389040e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7041245497695400
    v  = 0.4699877075860818e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6744033088306065e-1
    v  = 0.2099942281069176e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1678684485334166
    v  = 0.3172269150712804e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2793559049539613
    v  = 0.3832051358546523e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3935264218057639
    v  = 0.4252193818146985e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5052629268232558
    v  = 0.4513807963755000e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.6107905315437531
    v  = 0.4657797469114178e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1135081039843524
    b  = 0.3331954884662588e-1
    v  = 0.2733362800522836e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1612866626099378
    b  = 0.7247167465436538e-1
    v  = 0.3235485368463559e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2100786550168205
    b  = 0.1151539110849745
    v  = 0.3624908726013453e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2592282009459942
    b  = 0.1599491097143677
    v  = 0.3925540070712828e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3081740561320203
    b  = 0.2058699956028027
    v  = 0.4156129781116235e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3564289781578164
    b  = 0.2521624953502911
    v  = 0.4330644984623263e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4035587288240703
    b  = 0.2982090785797674
    v  = 0.4459677725921312e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4491671196373903
    b  = 0.3434762087235733
    v  = 0.4551593004456795e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4928854782917489
    b  = 0.3874831357203437
    v  = 0.4613341462749918e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5343646791958988
    b  = 0.4297814821746926
    v  = 0.4651019618269806e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5732683216530990
    b  = 0.4699402260943537
    v  = 0.4670249536100625e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2214131583218986
    b  = 0.3873602040643895e-1
    v  = 0.3549555576441708e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2741796504750071
    b  = 0.8089496256902013e-1
    v  = 0.3856108245249010e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3259797439149485
    b  = 0.1251732177620872
    v  = 0.4098622845756882e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3765441148826891
    b  = 0.1706260286403185
    v  = 0.4286328604268950e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4255773574530558
    b  = 0.2165115147300408
    v  = 0.4427802198993945e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4727795117058430
    b  = 0.2622089812225259
    v  = 0.4530473511488561e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5178546895819012
    b  = 0.3071721431296201
    v  = 0.4600805475703138e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5605141192097460
    b  = 0.3508998998801138
    v  = 0.4644599059958017e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6004763319352512
    b  = 0.3929160876166931
    v  = 0.4667274455712508e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3352842634946949
    b  = 0.4202563457288019e-1
    v  = 0.4069360518020356e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3891971629814670
    b  = 0.8614309758870850e-1
    v  = 0.4260442819919195e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4409875565542281
    b  = 0.1314500879380001
    v  = 0.4408678508029063e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4904893058592484
    b  = 0.1772189657383859
    v  = 0.4518748115548597e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5375056138769549
    b  = 0.2228277110050294
    v  = 0.4595564875375116e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5818255708669969
    b  = 0.2677179935014386
    v  = 0.4643988774315846e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6232334858144959
    b  = 0.3113675035544165
    v  = 0.4668827491646946e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4489485354492058
    b  = 0.4409162378368174e-1
    v  = 0.4400541823741973e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5015136875933150
    b  = 0.8939009917748489e-1
    v  = 0.4514512890193797e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5511300550512623
    b  = 0.1351806029383365
    v  = 0.4596198627347549e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5976720409858000
    b  = 0.1808370355053196
    v  = 0.4648659016801781e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6409956378989354
    b  = 0.2257852192301602
    v  = 0.4675502017157673e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5581222330827514
    b  = 0.4532173421637160e-1
    v  = 0.4598494476455523e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6074705984161695
    b  = 0.9117488031840314e-1
    v  = 0.4654916955152048e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6532272537379033
    b  = 0.1369294213140155
    v  = 0.4684709779505137e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6594761494500487
    b  = 0.4589901487275583e-1
    v  = 0.4691445539106986e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld2702!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 2702 point Lebedev angular grid.
"""
function ld2702!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.2998675149888161e-4
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.4077860529495355e-3
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.2065562538818703e-1
    v  = 0.1185349192520667e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5250918173022379e-1
    v  = 0.1913408643425751e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.8993480082038376e-1
    v  = 0.2452886577209897e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1306023924436019
    v  = 0.2862408183288702e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1732060388531418
    v  = 0.3178032258257357e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2168727084820249
    v  = 0.3422945667633690e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2609528309173586
    v  = 0.3612790520235922e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3049252927938952
    v  = 0.3758638229818521e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3483484138084404
    v  = 0.3868711798859953e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3908321549106406
    v  = 0.3949429933189938e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4320210071894814
    v  = 0.4006068107541156e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4715824795890053
    v  = 0.4043192149672723e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5091984794078453
    v  = 0.4064947495808078e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5445580145650803
    v  = 0.4075245619813152e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6072575796841768
    v  = 0.4076423540893566e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6339484505755803
    v  = 0.4074280862251555e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6570718257486958
    v  = 0.4074163756012244e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6762557330090709
    v  = 0.4077647795071246e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6911161696923790
    v  = 0.4084517552782530e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7012841911659961
    v  = 0.4092468459224052e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7064559272410020
    v  = 0.4097872687240906e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6123554989894765e-1
    v  = 0.1738986811745028e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1533070348312393
    v  = 0.2659616045280191e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2563902605244206
    v  = 0.3240596008171533e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3629346991663361
    v  = 0.3621195964432943e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.4683949968987538
    v  = 0.3868838330760539e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5694479240657952
    v  = 0.4018911532693111e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.6634465430993955
    v  = 0.4089929432983252e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1033958573552305
    b  = 0.3034544009063584e-1
    v  = 0.2279907527706409e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1473521412414395
    b  = 0.6618803044247135e-1
    v  = 0.2715205490578897e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1924552158705967
    b  = 0.1054431128987715
    v  = 0.3057917896703976e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2381094362890328
    b  = 0.1468263551238858
    v  = 0.3326913052452555e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2838121707936760
    b  = 0.1894486108187886
    v  = 0.3537334711890037e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3291323133373415
    b  = 0.2326374238761579
    v  = 0.3700567500783129e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3736896978741460
    b  = 0.2758485808485768
    v  = 0.3825245372589122e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4171406040760013
    b  = 0.3186179331996921
    v  = 0.3918125171518296e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4591677985256915
    b  = 0.3605329796303794
    v  = 0.3984720419937579e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4994733831718418
    b  = 0.4012147253586509
    v  = 0.4029746003338211e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5377731830445096
    b  = 0.4403050025570692
    v  = 0.4057428632156627e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5737917830001331
    b  = 0.4774565904277483
    v  = 0.4071719274114857e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2027323586271389
    b  = 0.3544122504976147e-1
    v  = 0.2990236950664119e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2516942375187273
    b  = 0.7418304388646328e-1
    v  = 0.3262951734212878e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3000227995257181
    b  = 0.1150502745727186
    v  = 0.3482634608242413e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3474806691046342
    b  = 0.1571963371209364
    v  = 0.3656596681700892e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3938103180359209
    b  = 0.1999631877247100
    v  = 0.3791740467794218e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4387519590455703
    b  = 0.2428073457846535
    v  = 0.3894034450156905e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4820503960077787
    b  = 0.2852575132906155
    v  = 0.3968600245508371e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5234573778475101
    b  = 0.3268884208674639
    v  = 0.4019931351420050e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5627318647235282
    b  = 0.3673033321675939
    v  = 0.4052108801278599e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5996390607156954
    b  = 0.4061211551830290
    v  = 0.4068978613940934e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3084780753791947
    b  = 0.3860125523100059e-1
    v  = 0.3454275351319704e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3589988275920223
    b  = 0.7928938987104867e-1
    v  = 0.3629963537007920e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4078628415881973
    b  = 0.1212614643030087
    v  = 0.3770187233889873e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4549287258889735
    b  = 0.1638770827382693
    v  = 0.3878608613694378e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5000278512957279
    b  = 0.2065965798260176
    v  = 0.3959065270221274e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5429785044928199
    b  = 0.2489436378852235
    v  = 0.4015286975463570e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5835939850491711
    b  = 0.2904811368946891
    v  = 0.4050866785614717e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6216870353444856
    b  = 0.3307941957666609
    v  = 0.4069320185051913e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4151104662709091
    b  = 0.4064829146052554e-1
    v  = 0.3760120964062763e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4649804275009218
    b  = 0.8258424547294755e-1
    v  = 0.3870969564418064e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5124695757009662
    b  = 0.1251841962027289
    v  = 0.3955287790534055e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5574711100606224
    b  = 0.1679107505976331
    v  = 0.4015361911302668e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5998597333287227
    b  = 0.2102805057358715
    v  = 0.4053836986719548e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6395007148516600
    b  = 0.2518418087774107
    v  = 0.4073578673299117e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5188456224746252
    b  = 0.4194321676077518e-1
    v  = 0.3954628379231406e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5664190707942778
    b  = 0.8457661551921499e-1
    v  = 0.4017645508847530e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6110464353283153
    b  = 0.1273652932519396
    v  = 0.4059030348651293e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6526430302051563
    b  = 0.1698173239076354
    v  = 0.4080565809484880e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6167551880377548
    b  = 0.4266398851548864e-1
    v  = 0.4063018753664651e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6607195418355383
    b  = 0.8551925814238349e-1
    v  = 0.4087191292799671e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld3074!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 3074 point Lebedev angular grid.
"""
function ld3074!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.2599095953754734e-4
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.3603134089687541e-3
    n += gen_oh!(2, a, b, v, n, x, y, z, w)
    v  = 0.3586067974412447e-3
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.1886108518723392e-1
    v  = 0.9831528474385880e-4
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4800217244625303e-1
    v  = 0.1605023107954450e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.8244922058397242e-1
    v  = 0.2072200131464099e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1200408362484023
    v  = 0.2431297618814187e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1595773530809965
    v  = 0.2711819064496707e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2002635973434064
    v  = 0.2932762038321116e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2415127590139982
    v  = 0.3107032514197368e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2828584158458477
    v  = 0.3243808058921213e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3239091015338138
    v  = 0.3349899091374030e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3643225097962194
    v  = 0.3430580688505218e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4037897083691802
    v  = 0.3490124109290343e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4420247515194127
    v  = 0.3532148948561955e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4787572538464938
    v  = 0.3559862669062833e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5137265251275234
    v  = 0.3576224317551411e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5466764056654611
    v  = 0.3584050533086076e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6054859420813535
    v  = 0.3584903581373224e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6308106701764562
    v  = 0.3582991879040586e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6530369230179584
    v  = 0.3582371187963125e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6718609524611158
    v  = 0.3584353631122350e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6869676499894013
    v  = 0.3589120166517785e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6980467077240748
    v  = 0.3595445704531601e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7048241721250522
    v  = 0.3600943557111074e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5591105222058232e-1
    v  = 0.1456447096742039e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1407384078513916
    v  = 0.2252370188283782e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2364035438976309
    v  = 0.2766135443474897e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3360602737818170
    v  = 0.3110729491500851e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.4356292630054665
    v  = 0.3342506712303391e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5321569415256174
    v  = 0.3491981834026860e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.6232956305040554
    v  = 0.3576003604348932e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.9469870086838469e-1
    b  = 0.2778748387309470e-1
    v  = 0.1921921305788564e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1353170300568141
    b  = 0.6076569878628364e-1
    v  = 0.2301458216495632e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1771679481726077
    b  = 0.9703072762711040e-1
    v  = 0.2604248549522893e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2197066664231751
    b  = 0.1354112458524762
    v  = 0.2845275425870697e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2624783557374927
    b  = 0.1750996479744100
    v  = 0.3036870897974840e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3050969521214442
    b  = 0.2154896907449802
    v  = 0.3188414832298066e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3472252637196021
    b  = 0.2560954625740152
    v  = 0.3307046414722089e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3885610219026360
    b  = 0.2965070050624096
    v  = 0.3398330969031360e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4288273776062765
    b  = 0.3363641488734497
    v  = 0.3466757899705373e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4677662471302948
    b  = 0.3753400029836788
    v  = 0.3516095923230054e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5051333589553359
    b  = 0.4131297522144286
    v  = 0.3549645184048486e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5406942145810492
    b  = 0.4494423776081795
    v  = 0.3570415969441392e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5742204122576457
    b  = 0.4839938958841502
    v  = 0.3581251798496118e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1865407027225188
    b  = 0.3259144851070796e-1
    v  = 0.2543491329913348e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2321186453689432
    b  = 0.6835679505297343e-1
    v  = 0.2786711051330776e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2773159142523882
    b  = 0.1062284864451989
    v  = 0.2985552361083679e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3219200192237254
    b  = 0.1454404409323047
    v  = 0.3145867929154039e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3657032593944029
    b  = 0.1854018282582510
    v  = 0.3273290662067609e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4084376778363622
    b  = 0.2256297412014750
    v  = 0.3372705511943501e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4499004945751427
    b  = 0.2657104425000896
    v  = 0.3448274437851510e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4898758141326335
    b  = 0.3052755487631557
    v  = 0.3503592783048583e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5281547442266309
    b  = 0.3439863920645423
    v  = 0.3541854792663162e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5645346989813992
    b  = 0.3815229456121914
    v  = 0.3565995517909428e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5988181252159848
    b  = 0.4175752420966734
    v  = 0.3578802078302898e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2850425424471603
    b  = 0.3562149509862536e-1
    v  = 0.2958644592860982e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3324619433027876
    b  = 0.7330318886871096e-1
    v  = 0.3119548129116835e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3785848333076282
    b  = 0.1123226296008472
    v  = 0.3250745225005984e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4232891028562115
    b  = 0.1521084193337708
    v  = 0.3355153415935208e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4664287050829722
    b  = 0.1921844459223610
    v  = 0.3435847568549328e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5078458493735726
    b  = 0.2321360989678303
    v  = 0.3495786831622488e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5473779816204180
    b  = 0.2715886486360520
    v  = 0.3537767805534621e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5848617133811376
    b  = 0.3101924707571355
    v  = 0.3564459815421428e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6201348281584888
    b  = 0.3476121052890973
    v  = 0.3578464061225468e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3852191185387871
    b  = 0.3763224880035108e-1
    v  = 0.3239748762836212e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4325025061073423
    b  = 0.7659581935637135e-1
    v  = 0.3345491784174287e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4778486229734490
    b  = 0.1163381306083900
    v  = 0.3429126177301782e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5211663693009000
    b  = 0.1563890598752899
    v  = 0.3492420343097421e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5623469504853703
    b  = 0.1963320810149200
    v  = 0.3537399050235257e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6012718188659246
    b  = 0.2357847407258738
    v  = 0.3566209152659172e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6378179206390117
    b  = 0.2743846121244060
    v  = 0.3581084321919782e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4836936460214534
    b  = 0.3895902610739024e-1
    v  = 0.3426522117591512e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5293792562683797
    b  = 0.7871246819312640e-1
    v  = 0.3491848770121379e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5726281253100033
    b  = 0.1187963808202981
    v  = 0.3539318235231476e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6133658776169068
    b  = 0.1587914708061787
    v  = 0.3570231438458694e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6515085491865307
    b  = 0.1983058575227646
    v  = 0.3586207335051714e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5778692716064976
    b  = 0.3977209689791542e-1
    v  = 0.3541196205164025e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6207904288086192
    b  = 0.7990157592981152e-1
    v  = 0.3574296911573953e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6608688171046802
    b  = 0.1199671308754309
    v  = 0.3591993279818963e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6656263089489130
    b  = 0.4015955957805969e-1
    v  = 0.3595855034661997e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    ld3470!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)

Compute the 3470 point Lebedev angular grid.
"""
function ld3470!(x::AbstractVector, y::AbstractVector, z::AbstractVector, w::AbstractVector)
    a = 0.0
    b = 0.0
    n = 0

    v  = 0.2040382730826330e-4
    n += gen_oh!(1, a, b, v, n, x, y, z, w)
    v  = 0.3178149703889544e-3
    n += gen_oh!(3, a, b, v, n, x, y, z, w)
    a  = 0.1721420832906233e-1
    v  = 0.8288115128076110e-4
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4408875374981770e-1
    v  = 0.1360883192522954e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7594680813878681e-1
    v  = 0.1766854454542662e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1108335359204799;
    v  = 0.2083153161230153e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1476517054388567
    v  = 0.2333279544657158e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.1856731870860615
    v  = 0.2532809539930247e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2243634099428821
    v  = 0.2692472184211158e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.2633006881662727
    v  = 0.2819949946811885e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3021340904916283
    v  = 0.2920953593973030e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3405594048030089
    v  = 0.2999889782948352e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.3783044434007372
    v  = 0.3060292120496902e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4151194767407910
    v  = 0.3105109167522192e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4507705766443257
    v  = 0.3136902387550312e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.4850346056573187
    v  = 0.3157984652454632e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5176950817792470
    v  = 0.3170516518425422e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5485384240820989
    v  = 0.3176568425633755e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6039117238943308
    v  = 0.3177198411207062e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6279956655573113
    v  = 0.3175519492394733e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6493636169568952
    v  = 0.3174654952634756e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6677644117704504
    v  = 0.3175676415467654e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6829368572115624
    v  = 0.3178923417835410e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.6946195818184121
    v  = 0.3183788287531909e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7025711542057026
    v  = 0.3188755151918807e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.7066004767140119
    v  = 0.3191916889313849e-3
    n += gen_oh!(4, a, b, v, n, x, y, z, w)
    a  = 0.5132537689946062e-1
    v  = 0.1231779611744508e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.1297994661331225
    v  = 0.1924661373839880e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.2188852049401307
    v  = 0.2380881867403424e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.3123174824903457
    v  = 0.2693100663037885e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.4064037620738195
    v  = 0.2908673382834366e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.4984958396944782
    v  = 0.3053914619381535e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.5864975046021365
    v  = 0.3143916684147777e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.6686711634580175
    v  = 0.3187042244055363e-3
    n += gen_oh!(5, a, b, v, n, x, y, z, w)
    a  = 0.8715738780835950e-1
    b  = 0.2557175233367578e-1
    v  = 0.1635219535869790e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1248383123134007
    b  = 0.5604823383376681e-1
    v  = 0.1968109917696070e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1638062693383378
    b  = 0.8968568601900765e-1
    v  = 0.2236754342249974e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2035586203373176
    b  = 0.1254086651976279
    v  = 0.2453186687017181e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2436798975293774
    b  = 0.1624780150162012
    v  = 0.2627551791580541e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2838207507773806
    b  = 0.2003422342683208
    v  = 0.2767654860152220e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3236787502217692
    b  = 0.2385628026255263
    v  = 0.2879467027765895e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3629849554840691
    b  = 0.2767731148783578
    v  = 0.2967639918918702e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4014948081992087
    b  = 0.3146542308245309
    v  = 0.3035900684660351e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4389818379260225
    b  = 0.3519196415895088
    v  = 0.3087338237298308e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4752331143674377
    b  = 0.3883050984023654
    v  = 0.3124608838860167e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5100457318374018
    b  = 0.4235613423908649
    v  = 0.3150084294226743e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5432238388954868
    b  = 0.4574484717196220
    v  = 0.3165958398598402e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5745758685072442
    b  = 0.4897311639255524
    v  = 0.3174320440957372e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.1723981437592809
    b  = 0.3010630597881105e-1
    v  = 0.2182188909812599e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2149553257844597
    b  = 0.6326031554204694e-1
    v  = 0.2399727933921445e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2573256081247422
    b  = 0.9848566980258631e-1
    v  = 0.2579796133514652e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2993163751238106
    b  = 0.1350835952384266
    v  = 0.2727114052623535e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3407238005148000
    b  = 0.1725184055442181
    v  = 0.2846327656281355e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3813454978483264
    b  = 0.2103559279730725
    v  = 0.2941491102051334e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4209848104423343
    b  = 0.2482278774554860
    v  = 0.3016049492136107e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4594519699996300
    b  = 0.2858099509982883
    v  = 0.3072949726175648e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4965640166185930
    b  = 0.3228075659915428
    v  = 0.3114768142886460e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5321441655571562
    b  = 0.3589459907204151
    v  = 0.3143823673666223e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5660208438582166
    b  = 0.3939630088864310
    v  = 0.3162269764661535e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5980264315964364
    b  = 0.4276029922949089
    v  = 0.3172164663759821e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.2644215852350733
    b  = 0.3300939429072552e-1
    v  = 0.2554575398967435e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3090113743443063
    b  = 0.6803887650078501e-1
    v  = 0.2701704069135677e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3525871079197808
    b  = 0.1044326136206709
    v  = 0.2823693413468940e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3950418005354029
    b  = 0.1416751597517679
    v  = 0.2922898463214289e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4362475663430163
    b  = 0.1793408610504821
    v  = 0.3001829062162428e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4760661812145854
    b  = 0.2170630750175722
    v  = 0.3062890864542953e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5143551042512103
    b  = 0.2545145157815807
    v  = 0.3108328279264746e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5509709026935597
    b  = 0.2913940101706601
    v  = 0.3140243146201245e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5857711030329428
    b  = 0.3274169910910705
    v  = 0.3160638030977130e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6186149917404392
    b  = 0.3623081329317265
    v  = 0.3171462882206275e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.3586894569557064
    b  = 0.3497354386450040e-1
    v  = 0.2812388416031796e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4035266610019441
    b  = 0.7129736739757095e-1
    v  = 0.2912137500288045e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4467775312332510
    b  = 0.1084758620193165
    v  = 0.2993241256502206e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4883638346608543
    b  = 0.1460915689241772
    v  = 0.3057101738983822e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5281908348434601
    b  = 0.1837790832369980
    v  = 0.3105319326251432e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5661542687149311
    b  = 0.2212075390874021
    v  = 0.3139565514428167e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6021450102031452
    b  = 0.2580682841160985
    v  = 0.3161543006806366e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6360520783610050
    b  = 0.2940656362094121
    v  = 0.3172985960613294e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4521611065087196
    b  = 0.3631055365867002e-1
    v  = 0.2989400336901431e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.4959365651560963
    b  = 0.7348318468484350e-1
    v  = 0.3054555883947677e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5376815804038283
    b  = 0.1111087643812648
    v  = 0.3104764960807702e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5773314480243768
    b  = 0.1488226085145408
    v  = 0.3141015825977616e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6148113245575056
    b  = 0.1862892274135151
    v  = 0.3164520621159896e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6500407462842380
    b  = 0.2231909701714456
    v  = 0.3176652305912204e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5425151448707213
    b  = 0.3718201306118944e-1
    v  = 0.3105097161023939e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.5841860556907931
    b  = 0.7483616335067346e-1
    v  = 0.3143014117890550e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6234632186851500
    b  = 0.1125990834266120
    v  = 0.3168172866287200e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6602934551848843
    b  = 0.1501303813157619
    v  = 0.3181401865570968e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6278573968375105
    b  = 0.3767559930245720e-1
    v  = 0.3170663659156037e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)
    a  = 0.6665611711264577
    b  = 0.7548443301360158e-1
    v  = 0.3185447944625510e-3
    n += gen_oh!(6, a, b, v, n, x, y, z, w)

    return nothing
end


"""
    isavailable(n::Integer)

Check if the Lebedev quadrature rule of order `n` is available.
"""
function isavailable(n::Integer)
    return haskey(rules, n)
end


"""
    availablerules()

Print a table showing the available Lebedev quadrature rules.
"""
function availablerules()
    println("order | points")
    println("------|-------")
    for rule in rules
        @printf("%5i | %6i\n",rule.first,rule.second)
    end
end


# function order2points(n::Integer)

#     # sanity checks
#     if iseven(n)
#         error("only odd-numbered orders are available")
#     end

#     if n < 3 || n > 131
#         error("only orders between 3 and 131 are available")
#     end

#     # determine the rule, ÷ does integer division
#     rule = (n-1)÷2

#     # this array maps the rule number to the number of grid points
#     points = [
#        6,   14,   26,   38,   50,   74,   86,  110,  146,  170,
#      194,  230,  266,  302,  350,  386,  434,  482,  530,  590,
#      650,  698,  770,  830,  890,  974, 1046, 1118, 1202, 1274,
#     1358, 1454, 1538, 1622, 1730, 1814, 1910, 2030, 2126, 2222,
#     2354, 2450, 2558, 2702, 2810, 2930, 3074, 3182, 3314, 3470,
#     3590, 3722, 3890, 4010, 4154, 4334, 4466, 4610, 4802, 4934,
#     5090, 5294, 5438, 5606, 5810]

#     return points[rule]
# end


"""
    lebedev_by_order(n::Integer)

Compute the Lebedev angular grid corresponding to order number `n`.
"""
function lebedev_by_order(n::Integer)
    if isavailable(n)
        return lebedev_by_points(rules[n])
    else
        throw(ArgumentError("Quadrature order not available"))
    end
end


"""
    lebedev_by_points(n::Integer)

Compute the `n` point Lebedev angular grid.
"""
function lebedev_by_points(n::Integer)

    # initialize arrays before checking `n` is valid just to avoid repeating
    # code inside the if block
    x = zeros(n)
    y = zeros(n)
    z = zeros(n)
    w = zeros(n)

    if n == 6
        ld0006!(x, y, z, w)
    elseif n == 14
        ld0014!(x, y, z, w)
    elseif n == 26
        ld0026!(x, y, z, w)
    elseif n == 38
        ld0038!(x, y, z, w)
    elseif n == 50
        ld0050!(x, y, z, w)
    elseif n == 74
        ld0074!(x, y, z, w)
    elseif n == 86
        ld0086!(x, y, z, w)
    elseif n == 110
        ld0110!(x, y, z, w)
    elseif n == 146
        ld0146!(x, y, z, w)
    elseif n == 170
        ld0170!(x, y, z, w)
    elseif n == 194
        ld0194!(x, y, z, w)
    elseif n == 230
        ld0230!(x, y, z, w)
    elseif n == 266
        ld0266!(x, y, z, w)
    elseif n == 302
        ld0302!(x, y, z, w)
    elseif n == 350
        ld0350!(x, y, z, w)
    elseif n == 434
        ld0434!(x, y, z, w)
    elseif n == 590
        ld0590!(x, y, z, w)
    elseif n == 770
        ld0770!(x, y, z, w)
    elseif n == 974
        ld0974!(x, y, z, w)
    elseif n == 1202
        ld1202!(x, y, z, w)
    elseif n == 1454
        ld1454!(x, y, z, w)
    elseif n == 1730
        ld1730!(x, y, z, w)
    elseif n == 2030
        ld2030!(x, y, z, w)
    elseif n == 2354
        ld2354!(x, y, z, w)
    elseif n == 2702
        ld2702!(x, y, z, w)
    elseif n == 3074
        ld3074!(x, y, z, w)
    elseif n == 3470
        ld3470!(x, y, z, w)
    else
        throw(ArgumentError("n is not a valid number of grid points"))
    end

    return x, y, z, w
end

