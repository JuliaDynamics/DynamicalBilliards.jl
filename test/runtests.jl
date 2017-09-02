using DynamicalBilliards
using Base.Test
# Test options:
printinfo = true
longtests = true

# Get tests options:
if haskey(ENV, "DYNAMICALBILLIARDS_PRINTTEST")
  if ENV["DYNAMICALBILLIARDS_PRINTTEST"] == "false"
    printinfo = false
  end
end

if haskey(ENV, "DYNAMICALBILLIARDS_LONGTEST")
  if ENV["DYNAMICALBILLIARDS_LONGTEST"] == "false"
    longtests = false
  end
end

if longtests
    const partnum = 500
else
    const partnum = 10
end

include("testfunctions.jl")
print("DynamicalBilliards tests started at: ")
print(Dates.format(now(), "HH:MM:s"), "\n")
t = time()

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
@test check_lyapunov_spectrum(partnum, printinfo = printinfo)

print("DynamicalBilliards tests ended (successfully) at: ")
println(Dates.format(now(), "HH:MM:s"))
t = time() - t
println("Total time required was:")
println(round(t, 3), " seconds, or ", round(t/60, 3), " minutes")
