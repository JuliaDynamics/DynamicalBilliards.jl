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
    const partnum = 1000
else
    const partnum = 10
end

include("testfunctions.jl")
print("DynamicalBilliards tests started at: ")
print(Dates.format(now(), "HH:MM:s"), "\n")

# Perform tests:
@test check_straight_sinai(partnum, printinfo = printinfo)
@test check_magnetic_sinai(partnum, printinfo = printinfo)
@test check_straight_sinai_periodic(partnum, printinfo = printinfo)
@test check_magnetic_sinai_periodic(partnum, printinfo = printinfo)
@test check_magnetic_pinned(partnum, printinfo = printinfo)
@test check_previous_obstacle(partnum, printinfo = printinfo)
@test check_raysplitting_omega(partnum, printinfo = printinfo)
@test check_raysplitting_periodic(partnum, printinfo = printinfo)
@test check_splitterwall(partnum, printinfo = printinfo)
@test check_random_sinai(partnum, printinfo = printinfo)
@test check_klein_magnetic(partnum, printinfo = printinfo)

print("DynamicalBilliards tests ended (successfully) at: ")
println(Dates.format(now(), "HH:MM:s"))
