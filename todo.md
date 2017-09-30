#TODO:
* make billiard table tuple, and compare (using benchmark and code_warntype, NOT profiling)
* add ffmpeg to plotting
* make the billiard plotting functions take advantage of `isperiodic` to result
  in a unified code?...

# ISSUES
* raysplit magnetic bigfloat is currently not supported.
* Raysplit mangetic sometimes freezes (never stops evolution).
* Increase performance. A standard evolve call has 100,000 allocations
  even though everything is done through SVectors..
