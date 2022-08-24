
module Lebedev

using Printf: @printf
using OrderedCollections: OrderedDict

# the implemented quadrature rules are stored as (order,points) pairs
const rules = OrderedDict(
     3  =>    6,
     5  =>   14,
     7  =>   26,
     9  =>   38,
    11  =>   50,
    13  =>   74,
    15  =>   86,
    17  =>  110,
    19  =>  146,
    21  =>  170,
    23  =>  194,
    25  =>  230,
    27  =>  266,
    29  =>  302,
    31  =>  350,
    35  =>  434,
    41  =>  590,
    47  =>  770,
    53  =>  974,
    59  => 1202,
    65  => 1454,
    71  => 1730,
    77  => 2030,
    83  => 2354,
    89  => 2702,
    95  => 3074,
    101 => 3470,
    107 => 3890,
    113 => 4334,
    119 => 4802,
    125 => 5294
)

include("sphere_lebedev_rule.jl")

export lebedev_by_points
export lebedev_by_order
export isavailable
export availablerules
export getavailablerules
export getavailableorders
export getavailablepoints

end
