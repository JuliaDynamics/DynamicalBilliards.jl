using DynamicalBilliards
using Base.Test
# Test options:
printinfo = false
longtests = false

# Get tests options:
if haskey(ENV, "DYNAMICALBILLIARDS_PRINTTEST")
  if ENV["DYNAMICALBILLIARDS_PRINTTEST"] == "true"
      printinfo = true
    end
end

if haskey(ENV, "DYNAMICALBILLIARDS_LONGTEST")
  if ENV["DYNAMICALBILLIARDS_LONGTEST"] == "true"
      longtests = true
    end
end

if longtests
    include("longtests.jl")
else
    include("shorttests.jl")
endst

# Perform tests:
@test check_straight_sinai(printinfo = printinfo)
@test check_magnetic_sinai(printinfo = printinfo)
@test check_straight_sinai_periodic(printinfo = printinfo)
@test check_magnetic_sinai_periodic(printinfo = printinfo)
@test check_magnetic_pinned(printinfo = printinfo)
@test check_previous_obstacle(printinfo = printinfo)
@test check_raysplitting_omega(printinfo = printinfo)
@test check_raysplitting_periodic(printinfo = printinfo)
@test check_splitterwall(printinfo = printinfo)
