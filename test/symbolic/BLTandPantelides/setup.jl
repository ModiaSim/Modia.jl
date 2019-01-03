# Desired:
#   import ModiaMath
#
# In order that ModiaMath need not to be defined in the user environment, it is included via Modia:
import Modia
import Modia.ModiaMath

include("../../../src/language/ModiaLogging.jl")
include("../../../src/language/Instantiation.jl")

include("../../../src/symbolic/BLTandPantelides/BLTandPantelidesUtilities.jl")
include("../../../src/symbolic/BLTandPantelides/BLTandPantelides.jl")
