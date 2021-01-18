
using Lebedev
using Test
using SpecialFunctions

# Essentially all functions generating the points and the weights (ldxxxx!(...))
# depend on gen_oh!(...), such that we do a thorough unit test on the latter.
# For the ldxxxx!(...) functions we do an integration test and do not mock the
# calls to gen_oh!(...), hence those tests depend on the output of gen_oh!(...).

@testset "gen_oh" begin
    # initialize these arrays to some length greater than really needed
    N = 122
    x = zeros(N)
    y = zeros(N)
    z = zeros(N)
    w = zeros(N)
    a = b = 0.0
    v = 0.5
    n = 0

    # code 1, write in x, y, z, and w from index n+1
    code = 1
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)

    # test x
    @test x[1] ==  1.0
    @test x[2] == -1.0
    @test all(x[ 3:6 ] .==  0.0 )
    # test y
    @test all( y[ [1:2 ; 5:6] ] .==  0.0 )
    @test y[3] ==  1.0
    @test y[4] == -1.0
    # test z
    @test all( z[ 1:4 ] .==  0.0 )
    @test z[5] ==  1.0
    @test z[6] == -1.0
    # code 1 sets 6 elements
    @test n == 6

    # code 2, write in x, y, z, and w from index n+1
    code = 2
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)
    ã = 1.0 / sqrt(2.0)
    
    # test x
    @test all( x[ 7:10 ] .==  0.0 )
    @test all( x[ [11:12 ; 15:16] ] .==  ã )
    @test all( x[ [13:14 ; 17:18] ] .== -ã )
    # test y
    @test all( y[ [7:8  ; 15 ; 17] ] .==  ã )
    @test all( y[ [9:10 ; 16 ; 18] ] .== -ã )
    @test all( y[ 11:14 ] .==  0.0 )
    # test z
    @test all( z[ [7 ;  9 ; 11 ; 13] ] .==  ã )
    @test all( z[ [8 ; 10 ; 12 ; 14] ] .== -ã )
    @test all( z[ 15:18 ] .==  0.0 )
    # code 2 sets 12 elements
    @test n == 18

    # code 3, write in x, y, z, and w from index n+1
    code = 3
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)
    ã = 1.0 / sqrt(3.0)

    # test x
    @test all( x[ 19:22 ] .==  ã )
    @test all( x[ 23:26 ] .== -ã )
    # test y
    @test all( y[ [19:20 ; 23:24] ] .==  ã )
    @test all( y[ [21:22 ; 25:26] ] .== -ã )
    # test z
    @test all( z[ [19;21 ; 23;25] ] .==  ã )
    @test all( z[ [20;22 ; 24;26] ] .== -ã )
    # code 3 sets 8 elements
    @test n == 26

    # code 4, write in x, y, z, and w from index n+1
    code = 4
    a = 0.5
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)
    b̃ = sqrt(1.0 - 2.0*a^2)

    # test x
    @test all( x[ [27:30 ; 35:38] ] .==  a )
    @test all( x[ [31:34 ; 39:42] ] .== -a )
    @test all( x[ [43 ; 45 ; 47 ; 49] ] .==  b̃ )
    @test all( x[ [44 ; 46 ; 48 ; 50] ] .== -b̃ )
    # test y
    @test all( y[ [27:28 ; 31:32 ; 43:46] ] .==  a )
    @test all( y[ [29:30 ; 33:34 ; 47:50] ] .== -a )
    @test all( y[ [35 ; 37 ; 39 ; 41] ] .==  b̃ )
    @test all( y[ [36 ; 38 ; 40 ; 42] ] .== -b̃ )
    # test z
    @test all( z[ [27 ; 29 ; 31 ; 33] ] .==  b̃ )
    @test all( z[ [28 ; 30 ; 32 ; 34] ] .== -b̃ )
    @test all( z[ [35:36 ; 39:40 ; 43:44 ; 47:48] ] .==  a )
    @test all( z[ [37:38 ; 41:42 ; 45:46 ; 49:50] ] .== -a )
    # code 4 sets 24 elements
    @test n == 50

    # code 5, write in x, y, z, and w from index n+1
    code = 5
    a = 0.5
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)
    b̃ = sqrt(1.0 - a^2)

    # test x
    @test all( x[ [51:52 ; 59:60] ] .==  a )
    @test all( x[ [53:54 ; 61:62] ] .== -a )
    @test all( x[ [55:56 ; 63:64] ] .==  b̃ )
    @test all( x[ [57:58 ; 65:66] ] .== -b̃ )
    @test all( x[ 67:74 ] .== 0.0 )
    # test y
    @test all( y[ [51 ; 53 ; 71:72] ] .==  b̃ )
    @test all( y[ [52 ; 54 ; 73:74] ] .== -b̃ )
    @test all( y[ 59:66 ] .== 0.0 )
    @test all( y[ [55 ; 57 ; 67:68] ] .==  a )
    @test all( y[ [56 ; 58 ; 69:70] ] .== -a )
    # test z
    @test all( z[ 51:58 ] .== 0.0 )
    @test all( z[ [59 ; 61 ; 67 ; 69] ] .==  b̃ )
    @test all( z[ [60 ; 62 ; 68 ; 70] ] .== -b̃ )
    @test all( z[ [63 ; 65 ; 71 ; 73] ] .==  a )
    @test all( z[ [64 ; 66 ; 72 ; 74] ] .== -a )
    # code 5 sets 24 elements
    @test n == 74

    # code 6, write in x, y, z, and w from index n+1
    code = 6
    a = 0.5
    b = 0.5
    n += Lebedev.gen_oh!(code, a, b, v, n, x, y, z, w)
    c̃ = sqrt(1.0 - a^2 - b^2)

    # test x
    @test all( x[ [ 75:78  ;  83:86 ] ] .==  a )
    @test all( x[ [ 79:82  ;  87:90 ] ] .== -a )
    @test all( x[ [ 91:94  ;  99:102] ] .==  b )
    @test all( x[ [ 95:98  ; 103:106] ] .== -b )
    @test all( x[ [107:110 ; 115:118] ] .==  c̃ )
    @test all( x[ [111:114 ; 119:122] ] .== -c̃ )
    # test y
    @test all( y[ [75:76 ; 79:80 ; 115:116 ; 119:120] ] .==  b )
    @test all( y[ [77:78 ; 81:82 ; 117:118 ; 121:122] ] .== -b )
    @test all( y[ [83:84 ; 87:88 ;  99:100 ; 103:104] ] .==  c̃ )
    @test all( y[ [85:86 ; 89:90 ; 101:102 ; 105:106] ] .== -c̃ )
    @test all( y[ [91:92 ; 95:96 ; 107:108 ; 111:112] ] .==  a )
    @test all( y[ [93:94 ; 97:98 ; 109:110 ; 113:114] ] .== -a )
    # test z
    @test all( z[ [75  ; 77  ; 79  ; 81  ; 91  ; 93  ; 95  ; 97] ] .==  c̃ )
    @test all( z[ [76  ; 78  ; 80  ; 82  ; 92  ; 94  ; 96  ; 98] ] .== -c̃ )
    @test all( z[ [83  ; 85  ; 87  ; 89  ; 107 ; 109 ; 111 ; 113] ] .==  b )
    @test all( z[ [84  ; 86  ; 88  ; 90  ; 108 ; 110 ; 112 ; 114] ] .== -b )
    @test all( z[ [99  ; 101 ; 103 ; 105 ; 115 ; 117 ; 119 ; 121] ] .==  a )
    @test all( z[ [100 ; 102 ; 104 ; 106 ; 116 ; 118 ; 120 ; 122] ] .== -a )
    # code 6 sets 48 elements
    @test n == 122
    
    # test that all weights are correctly set
    @test all(w .== v)
    
    # test invalid code
    code = 234
    @test_throws ArgumentError Lebedev.gen_oh!(code, a, b, v, 0, x, y, z, w)
end


# @testset "utilities" begin
    # # only odd orders available
    # @test_throws ErrorException Lebedev.order2points(4)
    # # min order is 3
    # @test_throws ErrorException Lebedev.order2points(1) 
    # # max order is 131
    # @test_throws ErrorException Lebedev.order2points(133) 
    
    # # test that each order gives the expected number of points
    # expected_points = [
    #    6,   14,   26,   38,   50,   74,   86,  110,  146,  170,
    #  194,  230,  266,  302,  350,  386,  434,  482,  530,  590, 
    #  650,  698,  770,  830,  890,  974, 1046, 1118, 1202, 1274,
    # 1358, 1454, 1538, 1622, 1730, 1814, 1910, 2030, 2126, 2222,
    # 2354, 2450, 2558, 2702, 2810, 2930, 3074, 3182, 3314, 3470,
    # 3590, 3722, 3890, 4010, 4154, 4334, 4466, 4610, 4802, 4934,
    # 5090, 5294, 5438, 5606, 5810]
    
    # for n = 3:2:131
    #     @test Lebedev.order2points(n) == expected_points[(n-1)÷2]
    # end
# end


# Exact formula for the integration of a polynomial over the surface of a sphere
# See G.Folland, The American Mathematical Monthly, Vol. 108, No. 5 (May, 2001), pp. 446-448.
function sphere(k::Integer, l::Integer, m::Integer)
    # if any exponent is odd, the integral is zero by symmetry
    if isodd(k) || isodd(l) || isodd(m)
        return 0.0
    # else calculate its value based on the beta function (made up by gamma functions)
    else
        bₖ = (k+1)/2
        bₗ = (l+1)/2
        bₘ = (m+1)/2
        return 2.0*gamma(bₖ)*gamma(bₗ)*gamma(bₘ)/gamma(bₖ+bₗ+bₘ)
    end
end


# test error handling of lebedev_by_points
@testset "lebedev error handling" begin
    @test_throws ArgumentError lebedev_by_order(100)
    @test_throws ArgumentError lebedev_by_points(31)
end


# I adapt the relative tolerance to the tightest possible such that all tests are passed
tol = 5e-14
last_order = 0
for order in keys(Lebedev.rules)
    @testset "lebedev order $order" begin
        for n = last_order:order
            for k = 0:n
                for l = 0:n-k
                    m = n - l - k
            
                    # calculate quadrature points
                    x,y,z,w = lebedev_by_order(order)
            
                    integral = 0.0
                    for i = 1:length(w)
                        integral += 4.0*pi * w[i] * x[i]^k * y[i]^l * z[i]^m
                    end
            
                    @test integral ≈ sphere(k,l,m) rtol = tol atol = tol
                end
            end
        end
    end
    global last_order = order
end

# include("ccall.jl")
