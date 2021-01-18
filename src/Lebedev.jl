
module Lebedev

using Printf: @printf
using OrderedCollections: OrderedDict

# Implemented quadrature rules available in the package, stored as (order,points) pairs
const rules = OrderedDict(
     3 =>   6,
     5 =>  14,
     7 =>  26,
     9 =>  38,
    11 =>  50,
    13 =>  74,
    15 =>  86,
    17 => 110,
    19 => 146,
    21 => 170,
    23 => 194,
    25 => 230,
    27 => 266,
    29 => 302,
    31 => 350,
    35 => 434,
    41 => 590,
    47 => 770,
    53 => 974
)

include("sphere_lebedev_rule.jl")

export lebedev_by_points
export lebedev_by_order
export isavailable
export availablerules

end
