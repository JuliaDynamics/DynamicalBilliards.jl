MAKE BILLIARD TABLE A TOUPLE.

Add to documentation: notice that even though all walls are finite, they are considered infinite when calculating the collision time. This means that the only billiard tables allowed for this package are **convex polygons** (link to wiki).

In raysplitting page add warning to reset billiard after use

#TODO:
* make billiard table tuple, and compare (using benchmark and code_warntype, NOT profiling)
* add ffmpeg to plotting
* lyapunovspectrum

# ISSUES
* Ray splitting is faster than normal (WTFFF????)
* Iterate over a Billiard Table in a type stable manner
* Index a billiard table in a type stable manner but with dynamic index
* hexagonal periodic billiard doesnt work
* raysplit magnetic is IMPOSSIBLY AND EXTREMELY slow sometimes. There are some
  problematic orbits for sure (because sometimes its fast).
  Or maybe there is some really slow regression happening during relocate.
  One really needs to debug this.
