# Run this in the REPL, not Juno!!!

using PkgBenchmark
bresult = benchmarkpkg("DynamicalBilliards"; retune = true)
