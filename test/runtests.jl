
using Lebedev
using Test


@testset "gen_oh" begin
    # initialize these arrays to some length greater than really needed
    N = 122
    x = zeros(N)
    y = zeros(N)
    z = zeros(N)
    w = zeros(N)
    a = b = 0.0
    v = 0.5

    # code 1, write in x, y, z, and w from index n+1
    code = 1
    n = 0
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)

    @test x[1] ==  1.0
    @test count(!iszero,x) == 2
    @test y[4] == -1.0
    @test count(!iszero,y) == 2
    @test z[5] ==  1.0
    @test count(!iszero,z) == 2
    # code 1 sets 6 elements
    @test n == 6

    # code 2, write in x, y, z, and w from index n+1
    code = 2
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)

    ã = 1.0 / sqrt(2.0)
    @test x[13] ≈ -ã
    @test count(!iszero,x[7:n]) == 8
    @test y[10] ≈ -ã
    @test count(!iszero,y[7:n]) == 8
    @test z[11] ≈  ã
    @test count(!iszero,z[7:n]) == 8
    # code 2 sets 12 elements
    @test n == 18

    # code 3, write in x, y, z, and w from index n+1
    code = 3
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)

    ã = 1.0 / sqrt(3.0)
    @test x[19] ≈  ã
    @test count(!iszero,x[19:n]) == 8
    @test y[20] ≈  ã
    @test count(!iszero,y[19:n]) == 8
    @test z[21] ≈  ã
    @test count(!iszero,z[19:n]) == 8
    # code 3 sets 8 elements
    @test n == 26

    # code 4, write in x, y, z, and w from index n+1
    code = 4
    a = 0.5
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)

    b̃ = sqrt(1.0 - 2.0*a^2)
    @test x[30] ≈  a
    @test x[50] ≈ -b̃
    @test count(!iszero,x[27:n]) == 24
    @test y[29] ≈ -a
    @test y[36] ≈ -b̃
    @test count(!iszero,y[27:n]) == 24
    @test z[27] ≈  b̃
    @test z[37] ≈ -a
    @test count(!iszero,z[27:n]) == 24
    # code 4 sets 24 elements
    @test n == 50

    # code 5, write in x, y, z, and w from index n+1
    code = 5
    a = 0.5
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)

    b̃ = sqrt(1.0 - a^2)
    @test x[51] ≈  a
    @test x[58] ≈ -b̃
    @test count(!iszero,x[51:n]) == 16
    @test y[56] ≈ -a
    @test y[73] ≈ -b̃
    @test count(!iszero,y[51:n]) == 16
    @test z[67] ≈  b̃
    @test z[72] ≈ -a
    @test count(!iszero,z[51:n]) == 16
    # code 5 sets 24 elements
    @test n == 74

    # code 6, write in x, y, z, and w from index n+1
    code = 6
    a = 0.5
    b = 0.5
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)

    c̃ = sqrt(1.0 - a^2 - b^2)
    @test x[75]  ≈  a
    @test x[104] ≈ -b
    @test x[114] ≈ -c̃
    @test count(!iszero,x[75:n]) == 48
    @test y[76]  ≈  b
    @test y[94]  ≈ -a
    @test y[105] ≈ -c̃
    @test count(!iszero,y[75:n]) == 48
    @test z[90]  ≈ -b
    @test z[97]  ≈  c̃
    @test z[122] ≈ -a
    @test count(!iszero,z[75:n]) == 48
    # code 6 sets 48 elements
    @test n == 122
    
    @test all(w .== v)
end
