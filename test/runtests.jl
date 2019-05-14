using DynamicalBilliards
using Test
using DynamicalBilliards.Testing

print("DynamicalBilliards tests start now.")
timetime = time()

include("basic_tests.jl")
include("extended_tests.jl")
include("various_tests.jl")
include("type_stability.jl")
include("phasespaces_tests.jl")
include("raysplt_tests.jl")
include("lyapunov.jl")

print("\nDynamicalBilliards tests ended (successfully)!")
timetime = time() - timetime
println("Total time required was:")
println(round(timetime, digits = 3), " seconds, or ", round(timetime/60, digits=3), " minutes")
