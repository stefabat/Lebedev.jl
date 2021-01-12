
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
        x[n+1] =   a; y[n+1] = 0.0; z[n+1] = 0.0; w[n+1] = v;
        x[n+2] =  -a; y[n+2] = 0.0; z[n+2] = 0.0; w[n+2] = v;
        x[n+3] = 0.0; y[n+3] =   a; z[n+3] = 0.0; w[n+3] = v;
        x[n+4] = 0.0; y[n+4] =  -a; z[n+4] = 0.0; w[n+4] = v;
        x[n+5] = 0.0; y[n+5] = 0.0; z[n+5] =   a; w[n+5] = v;
        x[n+6] = 0.0; y[n+6] = 0.0; z[n+6] =  -a; w[n+6] = v;
        n = 6
    elseif code == 2
        a = 1.0 / sqrt(2.0)
        x[n+ 1] = 0.0; y[n+ 1] =   a; z[n+ 1] =   a; w[n+ 1] = v;
        x[n+ 2] = 0.0; y[n+ 2] =   a; z[n+ 2] =  -a; w[n+ 2] = v;
        x[n+ 3] = 0.0; y[n+ 3] =  -a; z[n+ 3] =   a; w[n+ 3] = v;
        x[n+ 4] = 0.0; y[n+ 4] =  -a; z[n+ 4] =  -a; w[n+ 4] = v;
        x[n+ 5] =   a; y[n+ 5] = 0.0; z[n+ 5] =   a; w[n+ 5] = v;
        x[n+ 6] =   a; y[n+ 6] = 0.0; z[n+ 6] =  -a; w[n+ 6] = v;
        x[n+ 7] =  -a; y[n+ 7] = 0.0; z[n+ 7] =   a; w[n+ 7] = v;
        x[n+ 8] =  -a; y[n+ 8] = 0.0; z[n+ 8] =  -a; w[n+ 8] = v;
        x[n+ 9] =   a; y[n+ 9] =   a; z[n+ 9] = 0.0; w[n+ 9] = v;
        x[n+10] =   a; y[n+10] =  -a; z[n+10] = 0.0; w[n+10] = v;
        x[n+11] =  -a; y[n+11] =   a; z[n+11] = 0.0; w[n+11] = v;
        x[n+12] =  -a; y[n+12] =  -a; z[n+12] = 0.0; w[n+12] = v;
        n = 12
    elseif code == 3
        a = 1.0 / sqrt(3.0)
        x[n+1] =   a; y[n+1] =   a; z[n+1] =   a; w[n+1] = v;
        x[n+2] =   a; y[n+2] =   a; z[n+2] =  -a; w[n+2] = v;
        x[n+3] =   a; y[n+3] =  -a; z[n+3] =   a; w[n+3] = v;
        x[n+4] =   a; y[n+4] =  -a; z[n+4] =  -a; w[n+4] = v;
        x[n+5] =  -a; y[n+5] =   a; z[n+5] =   a; w[n+5] = v;
        x[n+6] =  -a; y[n+6] =   a; z[n+6] =  -a; w[n+6] = v;
        x[n+7] =  -a; y[n+7] =  -a; z[n+7] =   a; w[n+7] = v;
        x[n+8] =  -a; y[n+8] =  -a; z[n+8] =  -a; w[n+8] = v;
        n = 8
    elseif code == 4
        b = sqrt(1.0 - 2.0*a^2)
        x[n+ 1] =   a; y[n+ 1] =   a; z[n+ 1] =   b; w[n+ 1] = v;
        x[n+ 2] =   a; y[n+ 2] =   a; z[n+ 2] =  -b; w[n+ 2] = v;
        x[n+ 3] =   a; y[n+ 3] =  -a; z[n+ 3] =   b; w[n+ 3] = v;
        x[n+ 4] =   a; y[n+ 4] =  -a; z[n+ 4] =  -b; w[n+ 4] = v;
        x[n+ 5] =  -a; y[n+ 5] =   a; z[n+ 5] =   b; w[n+ 5] = v;
        x[n+ 6] =  -a; y[n+ 6] =   a; z[n+ 6] =  -b; w[n+ 6] = v;
        x[n+ 7] =  -a; y[n+ 7] =  -a; z[n+ 7] =   b; w[n+ 7] = v;
        x[n+ 8] =  -a; y[n+ 8] =  -a; z[n+ 8] =  -b; w[n+ 8] = v;
        x[n+ 9] =   a; y[n+ 9] =   b; z[n+ 9] =   a; w[n+ 9] = v;
        x[n+10] =   a; y[n+10] =  -b; z[n+10] =   a; w[n+10] = v;
        x[n+11] =   a; y[n+11] =   b; z[n+11] =  -a; w[n+11] = v;
        x[n+12] =   a; y[n+12] =  -b; z[n+12] =  -a; w[n+12] = v;
        x[n+13] =  -a; y[n+13] =   b; z[n+13] =   a; w[n+13] = v;
        x[n+14] =  -a; y[n+14] =  -b; z[n+14] =   a; w[n+14] = v;
        x[n+15] =  -a; y[n+15] =   b; z[n+15] =  -a; w[n+15] = v;
        x[n+16] =  -a; y[n+16] =  -b; z[n+16] =  -a; w[n+16] = v;
        x[n+17] =   b; y[n+17] =   a; z[n+17] =   a; w[n+17] = v;
        x[n+18] =  -b; y[n+18] =   a; z[n+18] =   a; w[n+18] = v;
        x[n+19] =   b; y[n+19] =   a; z[n+19] =  -a; w[n+19] = v;
        x[n+20] =  -b; y[n+20] =   a; z[n+20] =  -a; w[n+20] = v;
        x[n+21] =   b; y[n+21] =  -a; z[n+21] =   a; w[n+21] = v;
        x[n+22] =  -b; y[n+22] =  -a; z[n+22] =   a; w[n+22] = v;
        x[n+23] =   b; y[n+23] =  -a; z[n+23] =  -a; w[n+23] = v;
        x[n+24] =  -b; y[n+24] =  -a; z[n+24] =  -a; w[n+24] = v;
        n = 24
    elseif code == 5
        b = sqrt(1.0 - a^2)
        x[n+ 1] =   a; y[n+ 1] =   b; z[n+ 1] = 0.0; w[n+ 1] = v;
        x[n+ 2] =   a; y[n+ 2] =  -b; z[n+ 2] = 0.0; w[n+ 2] = v;
        x[n+ 3] =  -a; y[n+ 3] =   b; z[n+ 3] = 0.0; w[n+ 3] = v;
        x[n+ 4] =  -a; y[n+ 4] =  -b; z[n+ 4] = 0.0; w[n+ 4] = v;
        x[n+ 5] =   b; y[n+ 5] =   a; z[n+ 5] = 0.0; w[n+ 5] = v;
        x[n+ 6] =   b; y[n+ 6] =  -a; z[n+ 6] = 0.0; w[n+ 6] = v;
        x[n+ 7] =  -b; y[n+ 7] =   a; z[n+ 7] = 0.0; w[n+ 7] = v;
        x[n+ 8] =  -b; y[n+ 8] =  -a; z[n+ 8] = 0.0; w[n+ 8] = v;
        x[n+ 9] =   a; y[n+ 9] = 0.0; z[n+ 9] =   b; w[n+ 9] = v;
        x[n+10] =   a; y[n+10] = 0.0; z[n+10] =  -b; w[n+10] = v;
        x[n+11] =  -a; y[n+11] = 0.0; z[n+11] =   b; w[n+11] = v;
        x[n+12] =  -a; y[n+12] = 0.0; z[n+12] =  -b; w[n+12] = v;
        x[n+13] =   b; y[n+13] = 0.0; z[n+13] =   a; w[n+13] = v;
        x[n+14] =   b; y[n+14] = 0.0; z[n+14] =  -a; w[n+14] = v;
        x[n+15] =  -b; y[n+15] = 0.0; z[n+15] =   a; w[n+15] = v;
        x[n+16] =  -b; y[n+16] = 0.0; z[n+16] =  -a; w[n+16] = v;
        x[n+17] = 0.0; y[n+17] =   a; z[n+17] =   b; w[n+17] = v;
        x[n+18] = 0.0; y[n+18] =   a; z[n+18] =  -b; w[n+18] = v;
        x[n+19] = 0.0; y[n+19] =  -a; z[n+19] =   b; w[n+19] = v;
        x[n+20] = 0.0; y[n+20] =  -a; z[n+20] =  -b; w[n+20] = v;
        x[n+21] = 0.0; y[n+21] =   b; z[n+21] =   a; w[n+21] = v;
        x[n+22] = 0.0; y[n+22] =   b; z[n+22] =  -a; w[n+22] = v;
        x[n+23] = 0.0; y[n+23] =  -b; z[n+23] =   a; w[n+23] = v;
        x[n+24] = 0.0; y[n+24] =  -b; z[n+24] =  -a; w[n+24] = v;
        n = 24
    elseif code == 6
        c = sqrt(1.0 - a^2 - b^2)
        x[n+ 1] =   a; y[n+ 1] =   b; z[n+ 1] =   c; w[n+ 1] = v;
        x[n+ 2] =   a; y[n+ 2] =   b; z[n+ 2] =  -c; w[n+ 2] = v;
        x[n+ 3] =   a; y[n+ 3] =  -b; z[n+ 3] =   c; w[n+ 3] = v;
        x[n+ 4] =   a; y[n+ 4] =  -b; z[n+ 4] =  -c; w[n+ 4] = v;
        x[n+ 5] =  -a; y[n+ 5] =   b; z[n+ 5] =   c; w[n+ 5] = v;
        x[n+ 6] =  -a; y[n+ 6] =   b; z[n+ 6] =  -c; w[n+ 6] = v;
        x[n+ 7] =  -a; y[n+ 7] =  -b; z[n+ 7] =   c; w[n+ 7] = v;
        x[n+ 8] =  -a; y[n+ 8] =  -b; z[n+ 8] =  -c; w[n+ 8] = v;
        x[n+ 9] =   a; y[n+ 9] =   c; z[n+ 9] =   b; w[n+ 9] = v;
        x[n+10] =   a; y[n+10] =   c; z[n+10] =  -b; w[n+10] = v;
        x[n+11] =   a; y[n+11] =  -c; z[n+11] =   b; w[n+11] = v;
        x[n+12] =   a; y[n+12] =  -c; z[n+12] =  -b; w[n+12] = v;
        x[n+13] =  -a; y[n+13] =   c; z[n+13] =   b; w[n+13] = v;
        x[n+14] =  -a; y[n+14] =   c; z[n+14] =  -b; w[n+14] = v;
        x[n+15] =  -a; y[n+15] =  -c; z[n+15] =   b; w[n+15] = v;
        x[n+16] =  -a; y[n+16] =  -c; z[n+16] =  -b; w[n+16] = v;
        x[n+17] =   b; y[n+17] =   a; z[n+17] =   c; w[n+17] = v;
        x[n+18] =   b; y[n+18] =   a; z[n+18] =  -c; w[n+18] = v;
        x[n+19] =   b; y[n+19] =  -a; z[n+19] =   c; w[n+19] = v;
        x[n+20] =   b; y[n+20] =  -a; z[n+20] =  -c; w[n+20] = v;
        x[n+21] =  -b; y[n+21] =   a; z[n+21] =   c; w[n+21] = v;
        x[n+22] =  -b; y[n+22] =   a; z[n+22] =  -c; w[n+22] = v;
        x[n+23] =  -b; y[n+23] =  -a; z[n+23] =   c; w[n+23] = v;
        x[n+24] =  -b; y[n+24] =  -a; z[n+24] =  -c; w[n+24] = v;
        x[n+25] =   b; y[n+25] =   c; z[n+25] =   a; w[n+25] = v;
        x[n+26] =   b; y[n+26] =   c; z[n+26] =  -a; w[n+26] = v;
        x[n+27] =   b; y[n+27] =  -c; z[n+27] =   a; w[n+27] = v;
        x[n+28] =   b; y[n+28] =  -c; z[n+28] =  -a; w[n+28] = v;
        x[n+29] =  -b; y[n+29] =   c; z[n+29] =   a; w[n+29] = v;
        x[n+30] =  -b; y[n+30] =   c; z[n+30] =  -a; w[n+30] = v;
        x[n+31] =  -b; y[n+31] =  -c; z[n+31] =   a; w[n+31] = v;
        x[n+32] =  -b; y[n+32] =  -c; z[n+32] =  -a; w[n+32] = v;
        x[n+33] =   c; y[n+33] =   a; z[n+33] =   b; w[n+33] = v;
        x[n+34] =   c; y[n+34] =   a; z[n+34] =  -b; w[n+34] = v;
        x[n+35] =   c; y[n+35] =  -a; z[n+35] =   b; w[n+35] = v;
        x[n+36] =   c; y[n+36] =  -a; z[n+36] =  -b; w[n+36] = v;
        x[n+37] =  -c; y[n+37] =   a; z[n+37] =   b; w[n+37] = v;
        x[n+38] =  -c; y[n+38] =   a; z[n+38] =  -b; w[n+38] = v;
        x[n+39] =  -c; y[n+39] =  -a; z[n+39] =   b; w[n+39] = v;
        x[n+40] =  -c; y[n+40] =  -a; z[n+40] =  -b; w[n+40] = v;
        x[n+41] =   c; y[n+41] =   b; z[n+41] =   a; w[n+41] = v;
        x[n+42] =   c; y[n+42] =   b; z[n+42] =  -a; w[n+42] = v;
        x[n+43] =   c; y[n+43] =  -b; z[n+43] =   a; w[n+43] = v;
        x[n+44] =   c; y[n+44] =  -b; z[n+44] =  -a; w[n+44] = v;
        x[n+45] =  -c; y[n+45] =   b; z[n+45] =   a; w[n+45] = v;
        x[n+46] =  -c; y[n+46] =   b; z[n+46] =  -a; w[n+46] = v;
        x[n+47] =  -c; y[n+47] =  -b; z[n+47] =   a; w[n+47] = v;
        x[n+48] =  -c; y[n+48] =  -b; z[n+48] =  -a; w[n+48] = v;
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
    isavailable(n::Integer)

Check if the Lebedev quadrature rule of order `n` is available.
"""
function isavailable(n::Integer)
    return haskey(rules,n)
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
    elseif n == 302
        ld0302!(x, y, z, w)
    elseif n == 590
        ld0590!(x, y, z, w)
    else
        throw(ArgumentError("n is not a valid number of grid points"))
    end
    
    return x,y,z,w
end

