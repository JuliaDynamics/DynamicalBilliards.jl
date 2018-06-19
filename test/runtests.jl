using DynamicalBilliards
using Dates, LinearAlgebra
using Test
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
    elseif ENV["DYNAMICALBILLIARDS_LONGTEST"] == "true"
        longtests = true
    end
end

if longtests
    partnum = 1000
else
    partnum = 50
end

include("straight.jl")
include("magnetic.jl")
include("raysplit.jl")
include("various.jl")
include("lyapunov.jl")
include("psos.jl")

print("DynamicalBilliards tests started at: ")
print(Dates.format(now(), "HH:MM:s"), "\n")
t = time()

fnames = (
    straight_sinai, straight_periodic, magnetic_sinai, magnetic_periodic,
    type_stability, lyapunov_spectrum,
    lyapunov_magnetic, escape_times, stadium_bm, cut_psos, coordinates
    boundarymap_portion_test, phasespace_portion_test, meancoltimes
    raysplit_straight, raysplit_magnetic)

for f in fnames
    println()
    f(partnum, printinfo=printinfo)
end

print("\nDynamicalBilliards tests ended (successfully)")
t = time() - t
println("Total time required was:")
println(round(t, digits = 3), " seconds, or ", round(t/60, digits=3), " minutes")
