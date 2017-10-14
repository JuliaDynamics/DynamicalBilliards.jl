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
    elseif ENV["DYNAMICALBILLIARDS_LONGTEST"] == "true"
        longtests = true
    end
end

if longtests
    partnum = 500
else
    partnum = 10
end

include("straight.jl")
include("magnetic.jl")
include("raysplit.jl")
include("various.jl")
include("lyapunov.jl")

print("DynamicalBilliards tests started at: ")
print(Dates.format(now(), "HH:MM:s"), "\n")
t = time()

fnames = (
    straight_sinai, straight_periodic, magnetic_sinai, magnetic_periodic,
    raysplit_straight, raysplit_magnetic, type_stability, lyapunov_spectrum,
    escape_times)
for f in fnames
    println()
    f(partnum, printinfo=printinfo)
end

print("\nDynamicalBilliards tests ended (successfully) at: ")
println(Dates.format(now(), "HH:MM:s"))
t = time() - t
println("Total time required was:")
println(round(t, 3), " seconds, or ", round(t/60, 3), " minutes")
