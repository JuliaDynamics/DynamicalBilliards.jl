#TODO:
* make billiard table tuple, and compare (using benchmark and code_warntype, NOT profiling)
* add ffmpeg to plotting
* lyapunovspectrum
* make the billiard plotting functions take advantage of `isperiodic` to result
  in a unified code?...

# ISSUES
* Iterate over a Billiard Table in a type stable manner
* hexagonal periodic billiard doesnt work
* raysplit magnetic is IMPOSSIBLY AND EXTREMELY slow sometimes. There are some
  problematic orbits for sure (because sometimes its fast).
  Or maybe there is some really slow regression happening during relocate.
  One really needs to debug this.
* BigFloat and periodic walls does not work for magnetic propagation!
