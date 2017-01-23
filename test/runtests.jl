using DynamicalBilliards
using Base.Test
# Test options:
printinfo = false
longtests = false

# Get tests options:
if haskey(ENV, "DYNAMICALBILLIARDS_PRINTTEST")
  if ENV["DYNAMICALBILLIARDS_PRINTTEST") == true
      printinfo = true
    end
end

if haskey(ENV, "DYNAMICALBILLIARDS_LONGTEST")
  if ENV["DYNAMICALBILLIARDS_LONGTEST") == true
      printinfo = true
    end
end

if longtests
    include("longtests.jl")
else
    include("shorttests.jl")
end

# Perform tests:
@test check_straight_sinai_periodic(printinfo = printinfo)
