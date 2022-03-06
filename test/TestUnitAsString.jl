module TestUnitAsString

using Modia
using Test

v1 = 2.0u"m/s"
unit_v1 = unit(v1)
v1_unit = Modia.unitAsString( unit_v1  )   # = "m*s^-1"
v2_withoutUnit = 2.0
code = :( $v2_withoutUnit*@u_str($v1_unit) )  # = 2.0u"m*s^-1"
v2 = eval(code)
@show v1
@show v1_unit
@show v2

@test v1==v2
end