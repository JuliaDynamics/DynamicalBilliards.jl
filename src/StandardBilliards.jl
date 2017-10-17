using StaticArrays

export billiard_rectangle, billiard_sinai, billiard_polygon, billiard_lorentz,
billiard_raysplitting_showcase, billiard_hexagonal_sinai

const SV = SVector{2}
####################################################
## Famous/Standard Billiards
####################################################
"""
```julia
billiard_rectangle(x=1.0, y=1.0; setting = "standard")
```
Return a vector of obstacles that defines a rectangle billiard of size (`x`, `y`).

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type,
  enforcing periodicity at the boundaries
* "random" : The velocity is randomized upon collision.
* "ray-splitting" : All obstacles in the billiard table allow for ray-splitting.
"""
function billiard_rectangle(x=1.0, y=1.0; setting::String = "standard")

    x, y = promote(x,y)
    x = convert(AbstractFloat, x)
    bt = Obstacle{typeof(x)}[]
    o = typeof(x)(0.0)
    if setting == "standard"
        sp = [o,o]; ep = [o, y]; n = [x,o]
        leftw2 = InfiniteWall(sp, ep, n, "Left wall")
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw2 = InfiniteWall(sp, ep, n, "Right wall")
        sp = [o,y]; ep = [x, y]; n = [o,-y]
        topw2 = InfiniteWall(sp, ep, n, "Top wall")
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw2 = InfiniteWall(sp, ep, n, "Bottom wall")
        push!(bt, leftw2, rightw2, topw2, botw2)
    elseif setting == "periodic"
        sp = [o,o]; ep = [o, y]; n = [x,o]
        leftw = PeriodicWall(sp, ep, n, "Left periodic boundary")
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = PeriodicWall(sp, ep, n, "Right periodic boundary")
        sp = [o,y]; ep = [x, y]; n = [o,-y]
        topw = PeriodicWall(sp, ep, n, "Top periodic boundary")
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = PeriodicWall(sp, ep, n, "Bottom periodic boundary")
        push!(bt, leftw, rightw, topw, botw)
    elseif setting == "random"
        sp = [o,o]; ep = [o, y]; n = [x,o]
        leftw = RandomWall(sp, ep, n, "Left random wall")
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = RandomWall(sp, ep, n, "Right random wall")
        sp = [o,y]; ep = [x, y]; n = [o,-y]
        topw = RandomWall(sp, ep, n, "Top random wall")
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = RandomWall(sp, ep, n, "Bottom random wall")
        push!(bt, leftw, rightw, topw, botw)
    elseif setting == "ray-splitting"
        sp = [o,o]; ep = [o, y]; n = [x,o]
        leftw = SplitterWall(sp, ep, n, "Left ray-splitting wall")
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = SplitterWall(sp, ep, n, "Right ray-splitting wall")
        sp = [o,y]; ep = [x, y]; n = [o,-y]
        topw = SplitterWall(sp, ep, n, "Top ray-splitting wall")
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = SplitterWall(sp, ep, n, "Bottom ray-splitting wall")
        push!(bt, leftw, rightw, topw, botw)
    else
        throw(ArgumentError("The given setting=$setting is unknown."))
    end
    return bt
end

"""
```julia
billiard_sinai(r=0.25, x=1.0, y=1.0; setting = "standard")
```
Return a vector of obstacles that defines a Sinai billiard of size (`x`, `y`) with
a disk in its center, of radius `r`.

In the periodic case, the system is also known as "Lorentz Gas".

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type,
  enforcing periodicity at the boundaries
* "random" : The velocity is randomized upon collision.
* "ray-splitting" : All obstacles in the billiard table allow for ray-splitting.
"""
function billiard_sinai(r=0.25, x=1.0, y=1.0; setting = "standard")
    if (setting == "periodic") && (r>=x/2 || r>=y/2)
        es = "Disk radius too big for a periodic Sinai billiard.\n"
        es*= "Obstacles must not overlap with `PeriodicWall`s."
        error(es)
    end
    r, x, y = promote(r,x,y)
    bt = billiard_rectangle(x, y; setting = setting)
    c = [x/2, y/2]
    if setting == "random"
        centerdisk = RandomDisk(c, r, "Random disk")
    elseif setting == "ray-splitting"
        centerdisk = Antidot(c, r, "Antidot")
    else
        centerdisk = Disk(c, r, "Disk")
    end
    push!(bt, centerdisk)
end

"""
    billiard_lorentz(r=0.25, x=1.0, y=1.0)
Alias for `billiard_sinai(r,x,y; setting = "periodic")`.
"""
billiard_lorentz(r=0.25, x=1.0, y=1.0) = billiard_sinai(r,x,y; setting = "periodic")

"""
    billiard_polygon(n::Int, R, center = [0,0]; setting = "standard")
Return a vector of obstacles that defines a regular-polygonal billiard table
with `n` sides, radius `r` and given `center`.

Note: `R` denotes the so-called outer radius, not the inner one.

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type, enforcing periodicity
  at the boundaries. Only available for `n=4` or `n=6`.
* "random" : The velocity is randomized upon collision.
"""
function billiard_polygon(sides::Int, r::Real, center = [0,0]; setting = "standard")
    S = typeof(convert(AbstractFloat, r))
    bt = Obstacle{S}[]
    verteces = [S[r*cos(2π*i/sides), r*sin(2π*i/sides)] + center for i in 1:sides]

    if setting == "standard"
        T = InfiniteWall
        wallname = "wall"
    elseif setting == "periodic"
        if sides != 4 && sides != 6
            error("Polygonal and periodic billiard can exist only for `n=4` or `n=6`")
        end
        T = PeriodicWall
        wallname = "periodic wall"
        inr = S(r*cos(π/sides))
    elseif setting == "random"
        T = RandomWall
        wallname = "random wall"
    end

    for i in eachindex(verteces)
        N = length(verteces)
        starting = verteces[i]
        ending = verteces[mod1(i+1, N)]
        # Normal vector must look at where the particle is coming from
        w = ending - starting
        if setting == "periodic"
            normal = 2inr*normalize([-w[2], w[1]])
            wall = PeriodicWall(starting, ending, normal, wallname*" $i")
        else
            normal = [-w[2], w[1]]
            wall = InfiniteWall(starting, ending, normal, wallname*" $i")
        end
        push!(bt, wall)
    end
    return bt
end

"""
    billiard_hexagonal_sinai(r, R, center = [0,0]; setting = "standard")
Create a sinai-like billiard, which is a hexagon of outer radius `R`, containing
at its center (given by `center`) a disk of radius `r`. The `setting` keyword
is passed to `billiard_polygon`.
"""
function billiard_hexagonal_sinai(r::Real, R::Real, center = [0,0];
    setting = "standard")
    r, R = promote(r, R)
    T = typeof(r); center = T[center...]
    bt = billiard_polygon(6, R, center; setting = setting)
    DT = setting == "random" ? RandomDisk : Disk
    push!(bt, Disk(center, r))
    return bt
end




"""
    billiard_raysplitting_showcase(x=2.0, y=1.0, r1=0.3, r2=0.2) -> bt, rayspl
Showcase example billiard for ray-splitting processes. A rectangle `(x,y)` with a
SplitterWall at `x/2` and two disks at each side, with respective radii `r1`, `r2`.

**Notice**: This function returns a billiard table `bt` as well as a `rayspl`
dictionary!
"""
function billiard_raysplitting_showcase(x=2.0, y=1.0, r1=0.3, r2=0.2)

    r1≥x/4 || r2≥x/4 && throw(ArgumentError("Disks overlap with walls!"))

    sa = (θ, pflag, ω) -> pflag ? 2.0*θ : 0.5*θ
    Tp = (p) -> (θ, pflag, ω=0.0) -> begin
        if pflag
            abs(θ) < π/4 ? p*exp(-(θ)^2/2(π/8)^2) : 0.0
        else
            (1-p)*exp(-(θ)^2/2(π/4)^2)
        end
    end
    newo = ((x, bool) -> bool ? -0.5x : -2.0x)
    rayspl = Dict{Int, Vector{Function}}(
    5 => [Tp(0.5), sa, newo],
    6 => [Tp(0.35), sa, newo],
    7 => [Tp(0.65), sa, newo])

    bt = billiard_rectangle(x, y)
    sw = SplitterWall([x/2, 0.0], [x/2,y], [-1,0], true)
    push!(bt, sw)
    a1 = Antidot([x/4, y/2], r1, "Left Antidot")
    push!(bt, a1)
    a2 = Antidot([3x/4, y/2], r2, "Right Antidot")
    push!(bt, a2)

    return bt, rayspl
end

function billiard_square_mushroom(stem_length = 1.0, stem_width=0.2, cap_radious=1.0)

    sl = 1.0; sw = 0.2; cr =1.0
    leftcorn = SV(-sw/2, 0)
    rightcorn = SV(sw/2, 0)
    upleftcorn = SV(-sw/2, sl)
    uprightcorn = SV(sw/2, sl)

    S = typeof(convert(AbstractFloat, sl))
    bt = Obstacle{S}[]

    stembot = FiniteWall(leftcorn, rightcorn, SV(0, sw), true, "Stem bottom")
    stemleft = FiniteWall(leftcorn, upleftcorn, SV(sw, 0), false, "Stem left")
    stemright = FiniteWall(rightcorn, uprightcorn, SV(-sw, 0), false, "Stem right")

    push!(bt, stembot, stemleft, stemright)

    farleft = SV(-cr, sl)
    farright = SV(cr, sl)
    upfarleft = SV(-cr, sl+cr)
    upfarright = SV(cr, sl+cr)

    capbotleft = FiniteWall(upleftcorn, farleft, SV(0, sw), false)
    capleft = FiniteWall(farleft, upfarleft, SV(sw, 0), false)
    toptop = FiniteWall(upfarleft, upfarright, SV(0, -sw), false)
    capright = FiniteWall(farright, upfarright, SV(-sw, 0))
    capbotright = FiniteWall(uprightcorn, farright, SV(0, sw), false)

    push!(bt, capbotleft, capleft, toptop, capright, capbotright)

    return bt
end
