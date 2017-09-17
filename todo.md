MAKE BILLIARD TABLE A TOUPLE.

Add to documentation: notice that even though all walls are finite, they are considered infinite when calculating the collision time. This means that the only billiard tables allowed for this package are **convex polygons** (link to wiki).

In raysplitting page add warning to reset billiard after use

#TODO:
* rayplitting with magnetic
* make billiard table tuple, and compare (using benchmark and code_warntype, NOT profiling)
* restructure testset
* standardbilliards
* add ffmpeg to plotting
* lyapunovspectrum
* update docs LOL
* write some very general test for type stability,
  which will propagate with random sinai (to check these obstacles)

# ISSUES
* Ray splitting is faster than normal (WTFFF????)
* Iterate over a Billiard Table in a type stable manner
* Index a billiard table in a type stable manner but with dynamic index
